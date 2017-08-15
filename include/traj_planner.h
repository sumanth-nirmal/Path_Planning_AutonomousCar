# pragma once

#include <stdio.h>
#include <vector>
#include <eigen3/Eigen/Eigen>

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
  trajPlanner();
  ~trajPlanner();

  generateTrajctory(std::vector<double>& next_x_vals, std::vector<double>& next_y_vals);

private:
  laneNo getLane(double d);
  void minJerkTrajParam(double start[3], double end[3], double t, double &traj[6]);
};
