# Introduction

The repo is based on the particle filter implementation which I did in the TuBerlin 2019 Winter semester.

The problem description: given the provided map, we want to localize the robot.


## Procedure

1. uniformly distribute the particle set (pose and yaw)
2. Using odometry-based motion model to update
3. Implement likelihoodField Model

a. Precompute the distance map, distance to the nearest occupied grid

b. calculate the likelihood map based on the distance map

c. acquire the the probability where beam reachs at the likelihood map, accumulate the result of each beam. Based on that, we could update the weight.

4. Resample the particles

stochastic universal sampling

5. get the best hypothesis
