# Vehicle Path Follow Simulation in Ubuntu

## General description

* Create a path following simulation for vehicle using RK4, PID and Bicycle Model in C++
  * Minimalism code style, no 3rd parties libraries (except Eigen)

![](https://github.com/lenhatquang2512/Vehicle_Trajectory_Tracking_Sim/blob/main/gif/vehsim.gif)


User can choose any modes :

* Controller (P/PID/LQR)
* Discrete Propagation(Euler/RK4)
* Vehicle Dynamics (Naive/Advanced Bicycle model)
* Waypoint Generator (P2P/Sinusoidal/Cubic/Zigzag)
* Using/Not Using Low Pass Filter in PID
* Please modify which modes you want in the "config" Object
* Of course, for different scenario all gains need to be tuned again

## Requirement packages

* g++ in Linux **Ubuntu**

`sudo apt-get install g++`

* Gnuplot

`sudo apt-get install gnuplot`

* Eigen (**OPTIONAL**, *only installed if you want to work with LQR later*)

`sudo apt install libeigen3-dev`

## Usage

* Just compile normally as regular C++ code:

```sh
  g++  -std=c++14 -Wall -O3 -g VehiclePathFollowSim.cpp -o VehiclePathFollowSim
  ./VehiclePathFollowSim
```

* To change the controller type, in the **Config class** , modify:

```
  const bool usePID = true; //or false
```

* To change the waypoint trajectory, either toggle these 2 lines:

```
    const bool useZigZagWay = false; //or Sample/P2P
    const bool useSampleWay = true; //or Zigzag/P2P
```
* To change different propagation method (RK4/Euler) or Vehicle Dynamics Model(Naive/Advanced), change the enum :

```
  const PROPAGATOR_MODE propagator = RK4_NAIVE_DYNAMICS; 
  //Note it can be  RK4_NAIVE_DYNAMICS , EULER_NAIVE_DYNAMICS 
	// RK4_ADV_DYNAMICS , EULER_ADV_DYNAMICS 
```
* To use Low pass filter, change the value of beta (from 0 to 1, 0 means no use):

```
    const float beta = 0.9; // for low pass filter
```


## Notes

This is still under developing progress.
You can check this video that I performed : https://www.youtube.com/watch?v=G9c8e5KYjUQ


## References

* I use the Ricatti Solver in this repo : https://github.com/TakaHoribe/Riccati_Solver
* However, this LQR part is for future
