# pragma once

#include <stdio.h>
#include <vector>
#include <iostream>
#include <math.h>

// enum that describes the car location wrt to th ego car
enum carDir{
  LEFT_e = 0,
  RIGHT_e,
  FRONT_e,
  BEHIND_e

};

// enum for lanes
enum laneNo{
  NO_LANE_e = 0,
  LANE_1_e,
  LANE_2_e,
  LANE_3_e

};

class trajPlanner
{
public:

  void generateTrajctory(double car_x, double car_y, double car_yaw, std::vector<double>& next_x_vals, std::vector<double>& next_y_vals);
  trajPlanner();
  ~trajPlanner();

private:
  laneNo getLane(double d);
  double getDforLane(laneNo lane);
  void minJerkTrajParam(double start[3], double end[3], double t, std::vector<double> &traj);
};
