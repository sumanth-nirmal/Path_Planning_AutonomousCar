# pragma once

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <spline.h>
#include <chrono>
#include <unistd.h>

#define MetersPerSec_To_MPH          2.24
#define DESRIRED_VELOCITY_MPH        49.9
#define NO_OF_POINTS_PER_PATH        35
#define DISTANCE_THRESHOLD           20
#define HORIZON_X_THRESHOLD          20
#define DESIRED_ACCELERATION_MPS     10
#define MINIMUM_CAR_AHEAD_DISTANCE   3

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

// enum to acess the sensor fusion data
enum sensFusion{
  CAR_ID = 0,
  CAR_X,
  CAR_Y,
  CAR_VEL_X,
  CAR_VEL_Y,
  CAR_S,
  CAR_D
};

class trajPlanner
{
public:
  trajPlanner(std::string map_file);
  ~trajPlanner();

  void generateTrajctory(std::vector<double>& next_x_vals, std::vector<double>& next_y_vals);
  void update_ecar_params(double car_x, double car_y, double car_s, double car_d, double car_yaw, double car_speed);
  void update_previous_path(std::vector<double> previous_path_x, std::vector<double> previous_path_y, double end_path_s, double end_path_d);
  void update_sensor_fusion(std::vector<std::vector<double>> sensor_fusion);

private:

  // path for waypoints
  std::string wapypoints_path_;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  std::vector<double> map_waypoints_x_;
  std::vector<double> map_waypoints_y_;
  std::vector<double> map_waypoints_s_;
  std::vector<double> map_waypoints_dx_;
  std::vector<double> map_waypoints_dy_;

  // previous path
  std::vector<double> previous_path_x_;
  std::vector<double> previous_path_y_;
  double end_path_s_;
  double end_path_d_;

  // sensor data
  std::vector<std::vector<double>> sensor_fusion_;

  // ego car's data
  double car_x_;
  double car_y_;
  double car_s_;
  double car_d_;
  double car_yaw_;
  double car_speed_;

  // current lane
  laneNo curr_lane_;

  // current velocity
  double curr_vel_;

  double dt_;

  laneNo getLane(double d);
  double getDforLane(laneNo lane);
  void minJerkTrajParam(double start[3], double end[3], double t, std::vector<double> &traj);
  double distance(double x1, double y1, double x2, double y2);
  int ClosestWaypoint(double x, double y, std::vector<double> maps_x, std::vector<double> maps_y);
  int NextWaypoint(double x, double y, double theta, std::vector<double> maps_x, std::vector<double> maps_y);
  std::vector<double> getFrenet(double x, double y, double theta, std::vector<double> maps_x, std::vector<double> maps_y);
  std::vector<double> getXY(double s, double d, std::vector<double> maps_s, std::vector<double> maps_x, std::vector<double> maps_y);
  std::vector<laneNo> check_lanes(laneNo curr_lane);
  bool check_car_ahead(laneNo lane, double& car_ahead_vel, double &car_ahead_dist);
  void get_car_ahead(laneNo lane, double& car_ahead_vel, double &car_ahead_dist);
  bool lane_change_possible(laneNo lane, int fwd_th, int bwd_th);
  laneNo check_lane_change(void);
  laneNo get_intermediate_lane(laneNo dst_lane);
  void update_velocity(double desired_velocity);
};
