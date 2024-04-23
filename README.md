# Vehicle Path Follow Simulation in Ubuntu

## General description

* Create a path following simulation for vehicle using RK4, PID and Bicycle Model in C++
  * Minimalism code style, no 3rd parties libraries (except Eigen)

![](https://github.com/lenhatquang2512/Vehicle_Trajectory_Tracking_Sim/blob/main/gif/vehsim.gif)

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

## Notes

This is still under developing progress.

User can choose any modes :

* Controller (P/PID)
* Discrete Propagation(Euler/RK4)
* Vehicle Dynamics (Naive/Advanced Bicycle model)
* Waypoint Generator (P2P/Sinusoidal/Cubic/Zigzag)
* Please modify which modes you want in the "config" Object
* Of course, for different scenario all gains need to be tuned again

## References

* I use the Ricatti Solver in this repo : https://github.com/TakaHoribe/Riccati_Solver
* However, this LQR part is for future
