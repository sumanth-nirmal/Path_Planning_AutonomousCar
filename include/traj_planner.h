# pragma once

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <spline.h>


// For converting back and forth between radians and degrees.
inline constexpr double pi() { return M_PI; }
inline double deg2rad(double x) { return x * pi() / 180; }
inline double rad2deg(double x) { return x * 180 / pi(); }

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

  double distance(double x1, double y1, double x2, double y2);
  int ClosestWaypoint(double x, double y, std::vector<double> maps_x, std::vector<double> maps_y);
  int NextWaypoint(double x, double y, double theta, std::vector<double> maps_x, std::vector<double> maps_y);
  std::vector<double> getFrenet(double x, double y, double theta, std::vector<double> maps_x, std::vector<double> maps_y);
  std::vector<double> getXY(double s, double d, std::vector<double> maps_s, std::vector<double> maps_x, std::vector<double> maps_y);

  void generateTrajctory(double car_x, double car_y, double car_yaw, double car_s, double car_d, std::vector<double>& next_x_vals, std::vector<double>& next_y_vals);
  trajPlanner(std::string map_file);
  ~trajPlanner();

private:
  std::string wapypoints_path_;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  std::vector<double> map_waypoints_x_;
  std::vector<double> map_waypoints_y_;
  std::vector<double> map_waypoints_s_;
  std::vector<double> map_waypoints_dx_;
  std::vector<double> map_waypoints_dy_;

  laneNo getLane(double d);
  double getDforLane(laneNo lane);
  void minJerkTrajParam(double start[3], double end[3], double t, std::vector<double> &traj);
};
