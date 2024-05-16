/**
 * @file VehiclePathFollowSim.cpp
 * @author Quang Nhat Le (quangle@umich.edu)
 * @brief Complex trajectory tracking with advanced vehicle dynamics, propagated using RK4/Euler and Adaptive PID
 *        User can choose any modes : Controller (P/PID/LQR)
 *                                    Discrete Propagation(Euler/RK4)
 *                                    Vehicle Dynamics (Naive/Advanced Bicycle model)
 *                                    Waypoint Generator (P2P/Sinusoidal/Cubic/Zigzag)
 *        Please modify which modes you want in the "config" Object
 *        Of course, for different scenario all gains need to be tuned again
 *        For LQR , refer to : https://arxiv.org/html/2404.18312v1
 *                                
 * @version 3.0
 * @date 2024-05-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include <string>
#include <cstdio> // For the std::remove function
#include <memory>
#include <limits>
#include <random>
#include <time.h>

template<typename T>
using Waypoint = std::vector<std::array<T, 2>>;
template<typename T>
using Gain = std::array<T,2>;
template<typename T>
using Vsat = std::array<T,2>;

#define CLAMP(x, lower, upper) ((x) < (lower) ? (lower) : ((x) > (upper) ? (upper) : (x)))
#define PRINT_CMD(x) (std::cout << x << std::endl)
#define INF (std::numeric_limits<float>::infinity()) // T is float
#define PRINT_MAT(X) (std::cout << #X << ":\n" << X << std::endl << std::endl)

#define EIGEN_LIB 1

#ifdef EIGEN_LIB
    #include <Eigen/Dense>
#endif
typedef enum 
{
	RK4_NAIVE_DYNAMICS = 0,
	EULER_NAIVE_DYNAMICS = 1,
	RK4_ADV_DYNAMICS = 2,
	EULER_ADV_DYNAMICS = 3
}PROPAGATOR_MODE;

typedef enum{
    P_CONTROL = 0,
    PID_CONTROL = 1,
    LQR_CONTROL= 2
}CONTROLLER_ALG;

template<typename T>
class Config{
public:
    const unsigned int dimState = 3;
    const unsigned int dimControl = 2;
    const T dt = static_cast<T>(0.01);
    const std::string fpath = "vehicle_path.txt";
    const std::string fbody = "vehicle_body.txt";
    const std::string fwaypoint = "waypoints.txt";
    
    //-----------suggested tuned gains for Advanced Dynamics-------//
    // const Gain<T> Kp  = {{0.5,0.6}};
    // const Gain<T> Ki  = {{0.1,0.2}};
    // const Gain<T> Kd  = {{0.00,0.00}};
    // const Vsat<T> VClamp = {{-1.0,1.0}};
    // const Vsat<T> WClamp = {{-2,2}};

    //-----------suggested tuned gains for naive Dynamics-------//
    const Gain<T> Kp  = {{0.5,4.0}};
    const Gain<T> Ki  = {{0.2,1.0}};
    const Gain<T> Kd  = {{0.05,0.01}};
    const Vsat<T> VClamp = {{-3,3}};
    const Vsat<T> WClamp = {{-INF,INF}};

    const T goal_tol = 0.2;
    // const bool usePID = true; //or just P-control
    const bool useZigZagWay = false; //or Sample/P2P
    const bool useSampleWay = true; //or Zigzag/P2P
    const float beta = 0.1; // for low pass filter
    const PROPAGATOR_MODE propagator = RK4_NAIVE_DYNAMICS;
    const CONTROLLER_ALG controller = LQR_CONTROL;
    const bool LQRSaturated = false;
    // const float lqr_tol = 0.05;
    const unsigned int lqr_iter_max = 50;
};

#ifdef EIGEN_LIB
/* Itereation method for discrete model */
bool solveRiccatiIterationD(const Eigen::MatrixXf &Ad,
                            const Eigen::MatrixXf &Bd, const Eigen::MatrixXf &Q,
                            const Eigen::MatrixXf &R, Eigen::MatrixXf &P,
                            const float &tolerance = 1.E-5,
                            const unsigned int iter_max = 100000);
#endif

template<typename T>//from c++11 no need typedef,just warning
class STATE{ 
public:
    T x;
    T y;
    T yaw;

    STATE(void) = default;

    STATE(T x0, T y0, T yaw0):
        x(x0), y(y0), yaw(yaw0){}

    // Operator overload for addition of two STATEs
    STATE<T> operator+(const STATE<T>& rhs) const {
        return {x + rhs.x, y + rhs.y, yaw + rhs.yaw};
    }

    STATE<T> operator*(const T t) const {
        return {x * t, y * t, yaw * t};
    }

    // If modify the current STATE with the result of the addition
    STATE<T>& operator+=(const STATE<T>& rhs) {
        x += rhs.x;
        y += rhs.y;
        yaw += rhs.yaw;
        return *this;
    }

    #ifdef EIGEN_LIB
        //Overloading * operator for Eigen::Vector3f
        STATE<T> operator*(const Eigen::Vector3f &mat) const{
            return {x * mat(0), y * mat(1), yaw * mat(2)};
        }
    #endif

};
template<typename T>
class CONTROL{
public:
    T v;
    T w;
    CONTROL(void) = default;
};
//Function Pointer Prototypes
typedef STATE<float> (*DynamicsFunc)(const STATE<float>, const CONTROL<float>);

template<typename T>
constexpr T deg2rad(const T deg) {
    return static_cast<T>(deg * M_PI / 180.0);
}

template<typename T>
constexpr T rad2deg(const T rad) {
    return static_cast<T>(rad * 180.0/ M_PI);
}

template<typename T>
constexpr void nominalAngle(T &angle){ //between -pi and pi
    // angle = std::fmod(angle, 2.0f * M_PI);
    while (angle < -M_PI)
        angle += 2.0f * M_PI;
    while (angle > M_PI)
        angle -= 2.0f * M_PI;
}

template<typename T>
constexpr void clamp(T &u, const T umin, const T umax){
    if (u < umin) {
        u = umin;
    } else if (u > umax) {
        u = umax;
    }
}

template<typename T>
constexpr void rotate2D(T &xrot, T &yrot, const T x, const T y, const T angle){
	xrot = x * std::cos(angle) - y * std::sin(angle);
	yrot = x * std::sin(angle) + y * std::cos(angle);
}

template<typename T>
inline void HAL_Delay(const T sec){ //can not be constexpr
    usleep(sec * 1000000);
}

template<typename T>
inline T getRandomNum(const T randmin,const T randmax)
{
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<T> dis(randmin, randmax); 
//    std::uniform_int_distribution<T>(randmin,randmax);

   //different each invocation
   auto randNum = dis(gen);
   return randNum;
}

template<typename T>
constexpr T map(const T x, const T in_min, const T in_max, 
                    const T out_min, const T out_max){
    return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
}

// Helper function to compute the Euclidean distance between two points.
template<typename T>
constexpr T distance(const std::array<T, 2>& p1, const std::array<T, 2>& p2) {
    return std::hypot(p1[0] - p2[0], p1[1] - p2[1]);
}

template<typename T>
constexpr T digitalLowPassFilter(const T input, const T output,
                        const T coef){
    //this is Low pass exponential, coef aka beta is between 0 and 1
    return coef * input + (1-coef) * output;
}


template<typename T>
constexpr T magicPacejkaFormula(const T alpha, const T Fz, const T mu){
    const T B =  static_cast<T>(5.68); //from TireForce.csv, fitting
    const T C =  static_cast<T>(1.817);
    T Fy =  mu * Fz * std::sin(C * std::atan((B/mu) * alpha));
    return Fy;
}

template<typename T>
STATE<T> naiveDynamics(const STATE<T> X, const CONTROL<T> U){
    STATE<T> dX;
    dX.x = U.v * std::cos(X.yaw);
    dX.y = U.v * std::sin(X.yaw);
    dX.yaw = U.w;
    return dX;
}

template<typename T>
STATE<T> advancedDynamics(const STATE<T> X, const CONTROL<T> U){
    STATE<T> dX;
    T b = static_cast<T>(0.5f); // Distance from rear axle to vehicle's center of gravity
    T L = static_cast<T>(1.0f); // Wheelbase of the vehicle
    T phi = X.yaw;
    T u_in[2] = {U.v, U.w}; // Velocity and steering angle
    dX.x = u_in[0]*std::cos(phi) - u_in[0]* (b/L) * std::tan(u_in[1]) * std::sin(phi);
    dX.y = u_in[0]*std::sin(phi) + u_in[0]* (b/L) * std::tan(u_in[1]) * std::cos(phi); 
    dX.yaw = (1/L) *u_in[0] * std::tan(u_in[1]);
    return dX;
}

template<typename T>
STATE<T> forwardEuler(const STATE<T> X, const CONTROL<T> U, const T dt,DynamicsFunc dynamics){
    STATE<T> nextX ;
    STATE<T> dX = dynamics(X,U);
    // nextX.yaw =  X.yaw + dX.yaw * dt;
    // nextX.x = X.x + dX.x * dt; 
    // nextX.y = X.y + dX.y * dt;
    nextX = X;
    // nextX = X + dX * dt; //add using overloading + operator
    nextX += dX * dt;
    nominalAngle<T>(nextX.yaw);
    return nextX;
}

template<typename T>
STATE<T> RK4(const STATE<T> X, const CONTROL<T> U, const T dt, DynamicsFunc dynamics){
    STATE<T> nextX ;
    //Calculate the Runge-Kutta update for each state component
    STATE<T> k1 = dynamics(X,U);
    STATE<T> k2 = dynamics(X + k1 * (0.5*dt),U);
    STATE<T> k3 = dynamics(X + k2 * (0.5*dt),U);
    STATE<T> k4 = dynamics(X + k3 * (dt),U);
    //Update the state using the weighted sum of k1, k2, k3, and k4
    nextX = X + (k1+k2*2+k3*2+k4) * (dt/6);
    nominalAngle<T>(nextX.yaw);
    return nextX;
}

template<typename T>
void sampleCubicWaypointGenerator(Waypoint<T> &waypoints){
    for (T x = static_cast<T>(0.0); x <= static_cast<T>(3.0); x += static_cast<T>(0.5)) {
        // Calculate y using a cubic polynomial (e.g., y = ax^3 + bx^2 + cx+d)
        T a = 0.1;
        T b = 0.1;
        T c = 0.1;
        T d = 0.0;
        T y = a * std::pow(x,3) + b * std::pow(x,2)  + c * x + d;
        // Add the (x, y) point as a waypoint
        std::array<T,2> waypoint {x,y};
        waypoints.push_back(waypoint);
    }
}

template<typename T>
void sampleSineWaypointGenerator(Waypoint<T> &waypoints){
    T amplitude = 2.0;  // Amplitude of the sine wave
    T frequency = 1.0;  // How many waves in the 0 to 6 range
    T phase = 0.8;      // Phase shift, if any
    T yOffset = 1.0;    // Vertical offset to ensure the path doesn't go negative

    // create waypoints from x = 0 to 6 to match the range of the original zigzag
    for (T x = static_cast<T>(0.0); x <= static_cast<T>(6.0); x += static_cast<T>(0.5)) {
        T y = amplitude * std::sin(frequency * x + phase) + yOffset;
        waypoints.push_back({x, y});
    }
}

template<typename T>
void savePath(const std::string fname, const std::vector<STATE<T>> &path, const std::vector<T> &timeVec){
    std::ofstream MyFile(fname);

    for (size_t i = 0; i < path.size(); i++)
    {
        MyFile << path[i].x << " " << path[i].y << " " << path[i].yaw << " " << timeVec[i] << std::endl;
    }

    MyFile.close();
}

template<typename T>
void plotVehicle(std::ofstream &VehicleBody, const  STATE<T> X, const std::string fname){
    VehicleBody.open(fname, std::ofstream::out | std::ofstream::trunc);
    T l = static_cast<T>(0.5);
    T x1 = l/2; T y1 = l/2;
    T x2 = -l/2; T y2 = l/2;
    T x3 = -l/2; T y3 = -l/2;
    T x4 = l/2; T y4 = -l/2;

    T x11 = x1 * std::cos(X.yaw) - y1 * std::sin(X.yaw) + X.x;
    T y11 = x1 * std::sin(X.yaw) + y1 * std::cos(X.yaw) + X.y;
    T x22 = x2 * std::cos(X.yaw) - y2 * std::sin(X.yaw) + X.x;
    T y22 = x2 * std::sin(X.yaw) + y2 * std::cos(X.yaw) + X.y;
    T x33 = x3 * std::cos(X.yaw) - y3 * std::sin(X.yaw) + X.x;
    T y33 = x3 * std::sin(X.yaw) + y3 * std::cos(X.yaw) + X.y;
    T x44 = x4 * std::cos(X.yaw) - y4 * std::sin(X.yaw) + X.x;
    T y44 = x4 * std::sin(X.yaw) + y4 * std::cos(X.yaw) + X.y;

    VehicleBody << x11 << " " << y11 << std::endl;
    VehicleBody << x22 << " " << y22 << std::endl;
    VehicleBody << x33 << " " << y33 << std::endl;
    VehicleBody << x44 << " " << y44 << std::endl;
    VehicleBody << x11 << " " << y11 << std::endl;
    VehicleBody.close();

}

template<typename T>
void plotPath(std::ofstream &VehiclePath, const STATE<T> X, const CONTROL<T> U, const T time, const std::string fname){
    VehiclePath.open(fname, std::ofstream::out | std::ofstream::app); // Append mode
    VehiclePath << X.x << " " << X.y  << " " << X.yaw << " " << time  << U.v << " " << U.w << std::endl;
    VehiclePath.close(); // Close after writing
}

template<typename T>
void plotWaypoint(std::ofstream &WaypointFp, const Waypoint<T> waypoints, const std::string fname){
    WaypointFp.open(fname, std::ofstream::out | std::ofstream::trunc); // Truncate mode to overwrite
    for(const auto &waypoint : waypoints){
        WaypointFp << waypoint[0] << " " << waypoint[1] << std::endl;
    }
    WaypointFp.close(); // Close after writing
}

template<typename T>
void animationPlot(FILE *gp, const Config<T> config, const STATE<T> X){
    // Plot x versus y
    fprintf(gp, "set xrange [%f : %f]\n", X.x - 4.0f, X.x + 4.0f);
    fprintf(gp, "set yrange [%f : %f]\n", X.y - 4.0f, X.y + 4.0f);
    fprintf(gp, "set title 'Position (x vs y)'\n");
    // fprintf(gp, "plot \"%s\" using 1:2 with lines title 'x vs y'\n", (config.fpath).c_str());
    // fprintf(gp, "plot \"%s\" using 1:2 with lines title 'Path', \"%s\" using 1:2 with lines title 'Vehicle'\n",
    //        (config.fpath).c_str(), (config.fbody).c_str());
    // fprintf(gp, "plot \"%s\" using 1:2 with lines title 'Path', \"%s\" using 1:2 with filledcurves closed title 'Vehicle' fs solid 1.0 noborder\n",
    //             (config.fpath).c_str(), (config.fbody).c_str());  // to fill color
    fprintf(gp, "plot \"%s\" using 1:2 with lines title 'Path', \"%s\" using 1:2 with filledcurves closed title 'Vehicle' fs solid 1.0 noborder, \"%s\" using 1:2 with points pt 5 title 'Waypoints'\n",
            config.fpath.c_str(), config.fbody.c_str(), config.fwaypoint.c_str());
    fflush(gp);
}

template<typename T>
inline void PControl(T &U, const T error, const T Kp ){
    U = Kp * error;
}

template<typename T>
inline void PIDControl(T &U, const T error, const T prev_error,
       T &integral_error,const T Kp, const T Ki, const T Kd,
       const T dt,const T outMin, const T outMax,
       const STATE<T> prev_goal,
        const STATE<T> goal, const T beta, T &derivative_error){
    //Init integral error or any kind of error outside globally
    integral_error+= error * dt;
    clamp<T>(integral_error,outMin,outMax);
    // reset the integral if the reference is changed.
    if (goal.x != prev_goal.x || goal.y != prev_goal.y) {
        integral_error = static_cast<T>(0.0);
    }
    //  T derivative_error = (error - prev_error);
    T errorChange = error - prev_error;
    derivative_error = digitalLowPassFilter<T>(derivative_error,errorChange,beta);
    U = Kp * error + Ki * integral_error + Kd * (derivative_error/dt);
    clamp<T>(U,outMin,outMax);
     //need to update prev_error
}

template<typename T>
void computeError(const STATE<T> X, const STATE<T> goal, T &errV, T &errW){
    errV = static_cast<T>((std::pow(goal.x - X.x,2) + std::pow(goal.y - X.y,2)));
    errW = static_cast<T>(std::atan2(goal.y - X.y,goal.x - X.x));
    errW = errW - X.yaw;
    nominalAngle<T>(errW);
}

#ifdef EIGEN_LIB
bool solveRiccatiIterationD(const Eigen::MatrixXf &Ad,
                            const Eigen::MatrixXf &Bd, const Eigen::MatrixXf &Q,
                            const Eigen::MatrixXf &R, Eigen::MatrixXf &P,
                            const float &tolerance,
                            const unsigned int iter_max) {
  P = Q; // initialize

  Eigen::MatrixXf P_next;

  Eigen::MatrixXf AdT = Ad.transpose();
  Eigen::MatrixXf BdT = Bd.transpose();
  Eigen::MatrixXf Rinv = R.inverse();

  float diff;
  for (unsigned int i = 0; i < iter_max; ++i) {
    // -- discrete solver --
    P_next = AdT * P * Ad -
             AdT * P * Bd * (R + BdT * P * Bd).inverse() * BdT * P * Ad + Q;
    
    // PRINT_MAT(P_next);

    diff = (P_next - P).cwiseAbs().maxCoeff();
    P = P_next;
    // std::cout << "Iteration " << i << " with diff " << diff << std::endl;
    if (diff < tolerance) {
      std::cout << "iteration mumber = " << i << std::endl;
      return true;
    }
  }
  return false; // over iteration limit
}

template<typename T>
void naiveDiscreteZOHModel(const STATE<T> X, const CONTROL<T> U, const Config<T> config,
        Eigen::MatrixXf &Ad, Eigen::MatrixXf &Bd){
    Ad =  Eigen::MatrixXf::Zero(config.dimState, config.dimState);
    Bd =  Eigen::MatrixXf::Zero(config.dimState, config.dimControl);
    // Since the state does not change the system dynamics, the A matrix is 
    // the identity matrix -- no state feedback
    // Ad = I + A * dt
    Ad.setIdentity();
    // The B matrix for the discrete system
    // Position updates depend on velocity and yaw angle
    //Bd = B*dt
    Bd << std::cos(X.yaw) * config.dt, 0,
          std::sin(X.yaw) * config.dt, 0,
          0, config.dt;
}

template<typename T>
Eigen::VectorXf LQRFastControl(const STATE<T> errorState, 
                    const Eigen::MatrixXf Q, 
                    const Eigen::MatrixXf R, 
                    const Eigen::MatrixXf A, 
                    const Eigen::MatrixXf B, 
                    const Config<T> config) {
    // We want the system to stabilize at desired_state_xf.
    // Eigen::VectorXf x_error = actual_state_x - desired_state_xf;

    Eigen::MatrixXf AT = A.transpose();
    Eigen::MatrixXf BT = B.transpose();

    // Set the number of iterations
    const unsigned int N = config.lqr_iter_max;

    // Create an array of N + 1 Eigen::MatrixXf elements for P
    std::vector<Eigen::MatrixXf> P(N + 1);

    // Initialize to Q for the final cost
    P[N] = Q;

    // Compute the controller gains K.
    std::vector<Eigen::MatrixXf> K(N);

    // Dynamic programming to compute P
    for (unsigned int i = N; i > 0; --i) {
        // Discrete Algebraic Riccati equation
        P[i - 1] = Q + AT * P[i] * A - 
          (AT * P[i] * B) * (R + BT * P[i] * B).inverse() * (BT * P[i] * A);
    }

    // Compute K values
    for (unsigned int i = 0; i < N; ++i) {
        K[i] = (B.transpose() * P[i+1] * B + R).inverse() *B.transpose() * P[i+1] * A;
    }

    // You can also use the last computed K for the LQR input
    Eigen::VectorXf stateVec(config.dimState);
    stateVec << errorState.x, errorState.y, errorState.yaw;
    Eigen::VectorXf u_star = -K[N - 1] * stateVec;

    return u_star;
}


#endif

int main(int argc, char const *argv[])
{
    Config<float> config;
    // STATE X0 {.x = 0, .y = 0, .yaw = deg2rad(80) };
    // STATE goal {.x = 6, .y = 6, .yaw = deg2rad(0)};
    STATE<float> X0(0,0,deg2rad(0));
    STATE<float> goal(6,6,deg2rad(0));
    Waypoint<float> waypoints;
    if(config.useZigZagWay){
        waypoints = {{{0.0, 0.0}}, {{3.0, 0.0}}, {{3.0, 3.0}},{{6.0, 3.0}},{{6.0,6.0}}};
    }else if(config.useSampleWay){
        // sampleCubicWaypointGenerator(waypoints);
        sampleSineWaypointGenerator(waypoints);
    }else{
        waypoints = {{{6.0, 6.0}}}; // Point to Point (P2P) for easy Dynamic Debug
    }
    std::vector<float> timeVec;
    std::vector<STATE<float>> path;
    STATE<float> X {X0};
    CONTROL<float> U{.v=0,.w=0};
    bool Stop = false;
    float errV{0}; float errW{0};
    float integral_errV{0}; float integral_errW{0};
    float prev_errV{0}; float prev_errW{0};
    int count  = 0;
    FILE *gp;
    std::ofstream VehiclePath;
    std::ofstream VehicleBody;
    std::ofstream WaypointFp;
    goal.x = waypoints[0][0];
    goal.y = waypoints[0][1];
    size_t waypointSize = 0;
    STATE<float> prev_goal{goal}; //lastReference = reference;
    float derivative_error_V {0};
    float derivative_error_W {0};

    //LQR Setup
    #ifdef EIGEN_LIB
    Eigen::MatrixXf Q(3, 3); 
    Q.setIdentity();
    Q << 1.0 , 0 ,0,
        0 , 1.0, 0,
        0, 0 , 1.0;

    Eigen::MatrixXf R(2, 2);
    R.setIdentity();  // This is a placeholder. We can set actual R matrix values.
    R << 0.01, 0,
         0, 0.01;

    // Eigen::MatrixXf P = Eigen::MatrixXf::Zero(3, 3);  // Solution to Riccati equation
    
    // LQR Feedback gain matrix K
    // Eigen::MatrixXf K;
    Eigen::MatrixXf Ad;
    Eigen::MatrixXf Bd;
    //   STATE errorState(X0);

    #endif


    // Remove files if they exist
    std::remove(config.fpath.c_str()); // Deletes vehicle_path.txt
    std::remove(config.fbody.c_str()); // Deletes vehicle_body.txt
    std::remove(config.fwaypoint.c_str()); 
    
    // gp = popen("gnuplot","w");
    gp = popen("gnuplot -persist","w");

    if (gp == nullptr) {
        std::cerr << "Failed to open gnuplot." << std::endl;
        return 1;
    }else{
        fprintf(gp, "set colors classic\n");
        fprintf(gp, "set grid\n");
        fprintf(gp, "set xlabel x [m]\n");
        fprintf(gp, "set ylabel y [m]\n");
        fprintf(gp, "set tics font \"Arial, 14\"\n");
        fprintf(gp, "set xlabel font \"Arial, 14\"\n");
        fprintf(gp, "set ylabel font \"Arial, 14\"\n");
    }

    plotWaypoint<float>(WaypointFp,waypoints,config.fwaypoint);
    clock_t control_start_time = clock();

    while(!Stop){

        //------ Storing trajetory history----------------
        path.push_back(X);
        float time = count * config.dt;
        timeVec.push_back(time);
        
        //------------- Controller ---------------------
        computeError<float>(X,goal, errV,errW); //Error input for Controller block (aka setpoint - feedback)
        if(config.controller == PID_CONTROL){
            if(count==0) PRINT_CMD("USING FULL PID CONTROL");
            PIDControl<float>(U.v,errV,prev_errV,integral_errV,config.Kp[0],
             config.Ki[0],config.Kd[0],config.dt,config.VClamp[0],config.VClamp[1],
                        prev_goal,goal,config.beta,derivative_error_V);
            PIDControl<float>(U.w,errW,prev_errW,integral_errW,config.Kp[1],
             config.Ki[1],config.Kd[1],config.dt,-INF,INF,
             prev_goal,goal,config.beta,derivative_error_W);
            //Update pre_error
            prev_errV = errV;
            prev_errW = errW;
            //Update prev_goal
            prev_goal = goal;
        }else if(config.controller == P_CONTROL){
            if(count==0) PRINT_CMD("JUST USE P-CONTROL");
            PControl<float>(U.v,errV,config.Kp[0]);
            PControl<float>(U.w,errW,config.Kp[1]);
            // std::clamp(U.v,-1,1); //c++17
            clamp<float>(U.v,config.VClamp[0],config.VClamp[1]); // or U.v = CLAMP(U.v,-1,1);
            clamp<float>(U.w,config.WClamp[0],config.WClamp[1]);
        }else{
            #ifdef EIGEN_LIB
            //LQR
            if(count==0) PRINT_CMD("USING LQR CONTROL, Note that LQR is only used for NAIVE Dynamics");
            // Call function to get discrete linearized matrices (naive euler)
            naiveDiscreteZOHModel(X,U,config,Ad,Bd);
            // PRINT_CMD("Ad = ");
            // PRINT_MAT(Ad);
            PRINT_CMD("Bd = ");
            PRINT_MAT(Bd);
            STATE<float> errorState(X.x - goal.x, X.y - goal.y, errW);  // Construct an error state representation

            Eigen::VectorXf optimal_input = LQRFastControl(errorState, Q, R, Ad, Bd, config);
            U.v = optimal_input(0);
            U.w = optimal_input(1);

            // Saturate the inputs if necessary
            if(config.LQRSaturated){
                clamp(U.v, config.VClamp[0], config.VClamp[1]);
                clamp(U.w, config.WClamp[0],config.WClamp[1]);
            }
            #endif
        }


        //--------------Update motion of vehicle----------
        switch(config.propagator) {
            case (RK4_NAIVE_DYNAMICS):
                if(count==0) PRINT_CMD("RK4_NAIVE_DYNAMICS MODE");
                X = RK4<float>(X,U,config.dt,naiveDynamics);
                break;
            case (EULER_NAIVE_DYNAMICS):
                if(count==0) PRINT_CMD("EULER_NAIVE_DYNAMICS MODE");
                X = forwardEuler<float>(X,U,config.dt,naiveDynamics);
                break;
            case (RK4_ADV_DYNAMICS):
                if(count==0) PRINT_CMD("RK4_ADV_DYNAMICS MODE");
                X = RK4<float>(X,U,config.dt,advancedDynamics);
                break;
            case (EULER_ADV_DYNAMICS):
                if(count==0) PRINT_CMD("EULER_ADV_DYNAMICS MODE");
                X = forwardEuler<float>(X,U,config.dt,advancedDynamics);
                break;
            default:
                PRINT_CMD("Wrong propagator");
        }

        //Log terminal
        // std::cout << X.x << " " <<  X.y << " " << X.yaw << std::endl; // for Debug
        // std::cout << U.v << " " <<  U.w << std::endl;

        //-------Plot-----------------
        plotPath<float>(VehiclePath,X,U,time,config.fpath);
        plotVehicle<float>(VehicleBody,X,config.fbody);
        animationPlot<float>(gp,config,X);
        
        //------Terminate condition checking---------
        if(errV <= config.goal_tol){
            // U.v = 0; //No need for smooth transition
            // U.w = 0;
            waypointSize++;
            // Check if there are more waypoints
            if (waypointSize < waypoints.size()) {
                // Move to the next waypoint
                goal.x = waypoints[waypointSize][0];
                goal.y = waypoints[waypointSize][1];
            } else {
                // Stop the robot
                PRINT_CMD("Goal !");
                Stop=true;
                savePath<float>(config.fpath,path,timeVec);
                // break;
            }
        }
        count++;
        HAL_Delay<float>(config.dt);
    }

    clock_t control_end_time = clock();
    std::cout << "Computation time = " << static_cast<double>(control_end_time - control_start_time) / CLOCKS_PER_SEC
            << "sec." << std::endl;
    
    // VehiclePath.close();
    // VehicleBody.close();
    if(gp != NULL){
        fprintf(gp, "exit\n"); // Send an exit command to gnuplot
        pclose(gp);
    }

    // int retVal = system("killall -9 gnuplot\n");
    // int retVal = system("pkill -fx 'gnuplot -persist.*'");
    return 0;
}
