# Vehicle Path Follow Simulation in Ubuntu

## General description

* Create a path following simulation for vehicle using RK4, PID and Bicycle Model
  * Minimalism code style, no 3rd parties libraries (except Eigen)

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
