# Form-finding simulation of an elastic gridshell

# Overview

This project focused on the form-finding problem of an elastic gridshell.

<br/><img src='demo.gif' width="400">

To run this code, you should have a Linux Ubuntu system

# Make

g++ -fopenmp -I /usr/local/include/eigen-eigen-b3f3d4950030/ main.cpp world.cpp elasticRod.cpp elasticStretchingForce.cpp elasticBendingForce.cpp elasticTwistingForce.cpp externalGravityForce.cpp inertialForce.cpp contactForce.cpp timeStepper.cpp setInput.cpp -llapack -lGL -lglut -lGLU -Ofast -o gridShell

# Run 

./gridShell option.txt
