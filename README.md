# Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

---

## Problem Statement

Implement an Unscented Kalman filter to track a bicycle moving around the car using the Constant Turn Rate and Velocity (CTRV) motion model. Using radar and laser measurements develop a sensor fusion algorithm in C++ to track the bicycle's position,  velocity, yaw and yaw rate.

## Algorithm flow

The CTRV motion model is non-linear and the tracking is achieved using an Unscented Kalman filter. The main steps in programming a Kalman filter are:

1. **Initializing** Kalman filter variables
2. **Predicting** where the object is going to be after a time step Δt
3. **Updating** where the object is based on sensor measurements

The 2 measurements - laser and radar are used to track the bicycle using a Kalman filter and Unscented Kalman filter respectively.

## Output

The tracking accuracy is measured using the RMSE metric for position and velocity in x and y directions.

For obj_pose-laser-radar-synthetic-input.txt, the RMSE values are:

| Measurement | RMSE|
| --- | --- |
| Position (x) | 0.0610955 |
| Position (y) | 0.0858379 |
| Velocity (x) | 0.329842  |
| Velocity (y) | 0.213266  |

The position, velocity and yaw estimates vs the ground truth:

![Alt text](Pos_Est.png?raw=true)

![Alt text](Vel_Est.png?raw=true)

![Alt text](Yaw_Est.png?raw=true)

The Normalized Innovation Squared (NIS) for Laser and Radar measurements are captured:

![Alt text](Radar_NIS.png?raw=true)

![Alt text](Laser_NIS.png?raw=true)


## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it:  `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt ../output.txt`
