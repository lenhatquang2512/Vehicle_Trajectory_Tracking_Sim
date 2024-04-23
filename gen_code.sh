#!/bin/bash

clear

#compile
g++  -std=c++14 -Wall -O0 -g VehiclePathFollowSim.cpp -o VehiclePathFollowSim

#not necessary , make this file executable
chmod +x VehiclePathFollowSim

#run
./VehiclePathFollowSim