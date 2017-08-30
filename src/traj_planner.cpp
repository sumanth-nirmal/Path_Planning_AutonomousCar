#include <traj_planner.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"

using namespace std;

trajPlanner::trajPlanner(string map_file):
  wapypoints_path_(map_file)
{
  ifstream in_map_(wapypoints_path_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x_.push_back(x);
    map_waypoints_y_.push_back(y);
    map_waypoints_s_.push_back(s);
    map_waypoints_dx_.push_back(d_x);
    map_waypoints_dy_.push_back(d_y);
  }

  // start with middle lane
  curr_lane_ = LANE_2_e;

  // current velocity
  curr_vel_ = 0;

  // proc time
  dt_ = 0.02;
}

trajPlanner::~trajPlanner()
{

}

double trajPlanner::distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int trajPlanner::ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
  {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen)
    {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;
}

int trajPlanner::NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2( (map_y-y),(map_x-x) );

  double angle = abs(theta-heading);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  }

  return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> trajPlanner::getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> trajPlanner::getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}


laneNo trajPlanner::getLane(double d)
{
  laneNo ret = NO_LANE_e;

  if (d < 4)
  {
    ret = LANE_1_e;
  }
  else if (d >= 4 && d < 8)
  {
    ret = LANE_2_e;
  }
  else if (d >= 8)
  {
    ret = LANE_3_e;
  }

  return ret;
}

double trajPlanner::getDforLane(laneNo lane)
{
  double ret = 0;

  if (lane == LANE_1_e)
  {
    ret = 2;
  }
  else if (lane == LANE_2_e)
  {
    ret = 6;
  }
  else if (lane == LANE_3_e)
  {
    ret = 10;
  }

  return ret;
}

// source: http://mplab.ucsd.edu/tutorials/minimumJerk.pdf
// start and end is a vector with position, velocity and acceleration
void trajPlanner::minJerkTrajParam(double start[3], double end[3], double t, std::vector<double>& traj)
{
  // first 3 co-effiecients from the initial conditions
  double A0 = start[0];
  double A1 = start[1];
  double A2 = start[2]/2;

  // last 3 co - efficients from terminal conditions
  // minimum jerk in matrix form
  Eigen::Matrix3Xf T;

  T << std::pow(t, 3),      std::pow(t, 4),      std::pow(t, 5),
      3 * std::pow(t, 2),  4 * std::pow(t, 3),  5 * std::pow(t, 4),
      6 * t,               12 * std::pow(t, 2), 20 * std::pow(t, 3);

  Eigen::MatrixXf X(3,1);

  X << end[0] - A0 - A1*t - A2*t*t,
      end[1] - A1 - 2*A2*t,
      end[2] - 2*A2;

  Eigen::VectorXf A(3, 1);
  A = T.inverse() * X;

  // parameters for minimum jerk trajectory
  traj[0] = A0;   traj[1] = A1;     traj[2] = A2;
  traj[3] = A[0]; traj[4] = A[1];   traj[5] = A[2];
}

void trajPlanner::update_ecar_params(double car_x, double car_y, double car_s, double car_d, double car_yaw, double car_speed)
{
  car_x_ = car_x;
  car_y_ = car_y;
  car_s_ = car_s;
  car_d_ = car_d;
  car_yaw_ = car_yaw;
  car_speed_ = car_speed;
}

void trajPlanner::update_previous_path(vector<double> previous_path_x, vector<double> previous_path_y, double end_path_s, double end_path_d)
{
  previous_path_x_ = previous_path_x;
  previous_path_y_ = previous_path_y;
  end_path_s_ = end_path_s;
  end_path_d_ = end_path_d;
}

void trajPlanner::update_sensor_fusion(vector<vector<double>> sensor_fusion)
{
  sensor_fusion_.clear();

  for (int i=0; i<sensor_fusion.size(); i++)
  {
    sensor_fusion_.push_back(sensor_fusion[i]);
  }
}

vector<laneNo> trajPlanner::check_lanes(laneNo curr_lane)
{
  if(curr_lane == LANE_1_e)
  {
    return {LANE_2_e, LANE_3_e};
  }
  else if(curr_lane == LANE_2_e)
  {
    return {LANE_1_e, LANE_3_e};
  }
  else if(curr_lane == LANE_3_e)
  {
    return {LANE_2_e, LANE_1_e};
  }
}

// actually obselete
laneNo trajPlanner::get_intermediate_lane(laneNo dst_lane)
{
  if(getLane(car_d_) == LANE_1_e && dst_lane == LANE_3_e)
  {
    return LANE_2_e;
  }
  else if (getLane(car_d_) == LANE_3_e && dst_lane == LANE_1_e)
  {
    return LANE_2_e;
  }
}

bool trajPlanner::check_car_ahead(laneNo lane, double& car_ahead_vel, double& car_ahead_dist)
{
  bool ret = false;

  // loop thorugh the cars and check if they are nearer
  for (int i=0; i<sensor_fusion_.size(); i++)
  {
    // if car is in my lane
    if (lane == getLane(sensor_fusion_[i][CAR_D]))
    {
      // check if we can get too close to the car in the lane
      double sur_car_vel = std::sqrt(sensor_fusion_[i][CAR_VEL_X]*sensor_fusion_[i][CAR_VEL_X] + sensor_fusion_[i][CAR_VEL_Y]*sensor_fusion_[i][CAR_VEL_Y]);
      double sur_car_s = sensor_fusion_[i][CAR_S] + (previous_path_x_.size() * 0.02 * sur_car_vel);

      if ((sur_car_s > car_s_) &&    // if s neighbouring car is greater
          (sur_car_s - car_s_ < DISTANCE_THRESHOLD))
      {
        ret = true;
        car_ahead_vel = sur_car_vel;
        car_ahead_dist = sur_car_s - car_s_;
      }
    }
  }

  return ret;
}

void trajPlanner::get_car_ahead(laneNo lane, double& car_ahead_vel, double& car_ahead_dist)
{
  // initial value
  double low_car_s = std::numeric_limits<double>::max();

  // loop thorugh the cars and check if they are nearer
  for (int i=0; i<sensor_fusion_.size(); i++)
  {
    // if car is in my lane
    if (lane == getLane(sensor_fusion_[i][CAR_D]))
    {
      double sur_car_vel = std::sqrt(sensor_fusion_[i][CAR_VEL_X]*sensor_fusion_[i][CAR_VEL_X] + sensor_fusion_[i][CAR_VEL_Y]*sensor_fusion_[i][CAR_VEL_Y]);
      double sur_car_s = sensor_fusion_[i][CAR_S] + (previous_path_x_.size() * 0.02 * sur_car_vel);

      if (sur_car_s > car_s_ &&
          sur_car_s < low_car_s)
      {
        low_car_s = sur_car_s;

        // update the ahead velocity and the distance
        car_ahead_vel = sur_car_vel;
        car_ahead_dist = low_car_s - car_s_;
      }
    }
  }
}

bool trajPlanner::lane_change_possible(laneNo lane, int fwd_th, int bwd_th)
{
  bool ret = false;
  bool is_nocar_in_cutrange = true;

  std::cout << "checking in lane: " << lane << "\n";

  // loop thorugh the cars and check if they are nearer
  for (int i=0; i<sensor_fusion_.size(); i++)
  {
    // if car is in my lane
    if (lane == getLane(sensor_fusion_[i][CAR_D]))
    {
      // check if we can get too close to the car in the lane
      double sur_car_vel = std::sqrt(sensor_fusion_[i][CAR_VEL_X]*sensor_fusion_[i][CAR_VEL_X] + sensor_fusion_[i][CAR_VEL_Y]*sensor_fusion_[i][CAR_VEL_Y]);
      double sur_car_s = sensor_fusion_[i][CAR_S] + (previous_path_x_.size() * 0.02 * sur_car_vel);

      std::cout << "our car_s: " << car_s_ << " sur_car_s: " << sur_car_s << "\n";
      if (((sur_car_s > car_s_) && (sur_car_s - car_s_ < fwd_th)) ||
          ((sur_car_s < car_s_) && (car_s_ - sur_car_s < bwd_th)))
      {
        is_nocar_in_cutrange = false;
      }
    }
  }

  if (is_nocar_in_cutrange)
  {
    ret = true;
  }

  return ret;
}

laneNo trajPlanner::check_lane_change(void)
{
  laneNo ret_lane = NO_LANE_e;
  bool is_lane_op1_possible = false, is_lane_op2_possible = false;

  // get the lanes possible for opening
  std::vector<laneNo> lane = check_lanes(getLane(car_d_));

  if (lane_change_possible(lane[0], 20, 15))
  {
    ret_lane = lane[0];
    is_lane_op1_possible = true;
  }
  if(lane_change_possible(lane[1], 20, 15))
  {
    ret_lane = lane[1];
    is_lane_op2_possible = true;
  }

  // when lane change is possible in multiple lanes, change to the fastest lane
  if (is_lane_op1_possible == true &&
      is_lane_op2_possible == true)
  {
    double car_dist_lane_op1, car_dist_lane_op2, car_vel_lane_op1, car_vel_lane_op2;
    get_car_ahead(lane[0], car_vel_lane_op1, car_dist_lane_op1);
    get_car_ahead(lane[1], car_vel_lane_op2, car_dist_lane_op2);

    std::cout << "multiple lane change possible ";

    if (car_dist_lane_op1 > car_dist_lane_op2 &&
        car_vel_lane_op1 > car_vel_lane_op2)
    {
      ret_lane = lane[0];
      std::cout << "1\n";
    }
    else if (car_dist_lane_op2 > car_dist_lane_op1 &&
             car_vel_lane_op2 > car_vel_lane_op1)
    {
      ret_lane = lane[1];
      std::cout << "2\n";
    }
    else if (car_vel_lane_op2 > car_vel_lane_op1)
    {
      ret_lane = lane[1];
      std::cout << "3\n";
    }
    else if (car_vel_lane_op1 > car_vel_lane_op2)
    {
      ret_lane = lane[0];
      std::cout << "4\n";
    }
    else if (car_dist_lane_op2 > car_dist_lane_op1)
    {
      ret_lane = lane[1];
      std::cout << "5\n";
    }
    else if (car_dist_lane_op1 > car_dist_lane_op2)
    {
      ret_lane = lane[0];
      std::cout << "6\n";
    }
    else
    {
      ret_lane = lane[0];
      std::cout << "fall through case\n";
    }
  }

  return ret_lane;
}

void trajPlanner::update_velocity(double desired_velocity)
{
  if (car_speed_ < desired_velocity)
  {
    curr_vel_ += 1.3; //10*0.02*MetersPerSec_To_MPH;

    if (curr_vel_ > DESRIRED_VELOCITY_MPH)
    {
      curr_vel_ = DESRIRED_VELOCITY_MPH;
    }

    //std::cout << "curr vel: " << curr_vel_ << " car_speed: " << car_speed_ << "\n";
  }
  else
  {
    curr_vel_ -= 1.6; //10*0.02*MetersPerSec_To_MPH;
  }
}

void trajPlanner::generateTrajctory(std::vector<double>& next_x_vals, std::vector<double>& next_y_vals)
{
  vector<double> ptx, pty;
  double desire_vel;

  //reference car pose
  double ref_x = car_x_;
  double ref_y = car_y_;
  double ref_yaw = deg2rad(car_yaw_);

  // if there are previous points
  if (previous_path_x_.size() > 0)
  {
    // update the current car s position to the end point position of previous path
    car_s_ = end_path_s_;
  }

  // check if there is a car infornt which is too close
  double car_ahead_vel = 0, car_ahead_dist = MINIMUM_CAR_AHEAD_DISTANCE;  // default values if car ahead is not in the threshold
  if (check_car_ahead(getLane(car_d_), car_ahead_vel, car_ahead_dist))
  {
    std::cout << "\ncar ahead\n";
    laneNo lane = check_lane_change();
    std::cout << "shift to " << lane << "\n";

    // if sufficient place is there in front for lane change
    if (car_ahead_dist >= MINIMUM_CAR_AHEAD_DISTANCE)
    {
      if (lane != NO_LANE_e && abs(lane-getLane(car_d_)) == 1)
      {
        curr_lane_ = lane;
        desire_vel = DESRIRED_VELOCITY_MPH;
        std::cout << "shifted to " << curr_lane_ << "\n";
      }
      else if (lane != NO_LANE_e && abs(lane-getLane(car_d_)) > 1)
      {
        // check if there is possibility for mutiple lane switches
        if(lane_change_possible(get_intermediate_lane(lane), 15, 15))
        {
          curr_lane_ = LANE_2_e;
          desire_vel = DESRIRED_VELOCITY_MPH;
          std::cout << "inter lane shifted to " << curr_lane_ << "\n";
        }
        else
        {
          //set desired speed to the car in front
          desire_vel = car_ahead_vel;
          std::cout << "lane change not feasible" << lane << " car vel: " << car_ahead_vel << "\n";
        }
      }
      else
      {
        //set desired speed to the car in front
        desire_vel = car_ahead_vel;
        std::cout << "lane change not feasible" << lane << " car vel: " << car_ahead_vel << "\n";
      }
    }
  }
  else
  {
    desire_vel = DESRIRED_VELOCITY_MPH;
  }

  // update velocity
  update_velocity(desire_vel);

  // emergency stopping behaviour
  if (car_ahead_dist < MINIMUM_CAR_AHEAD_DISTANCE)
  {
    std::cout << "!!!!!..emergency breaking..!!!!!!\n";
    //curr_vel_ = 0;
    desire_vel = 0;
  }

  // if previous points are less than 2, (thats almost empty)
  if (previous_path_x_.size() < 2)
  {
    //we use the car state
    ptx.push_back(car_x_ - cos(car_yaw_));
    pty.push_back(car_y_ - sin(car_yaw_));

    ptx.push_back(car_x_);
    pty.push_back(car_y_);
  }
  // if it has points
  else
  {
    // we use the last 2 points from the previous points
    // also update reference pose
    ref_x = previous_path_x_[previous_path_x_.size() - 1];
    ref_y = previous_path_y_[previous_path_y_.size() - 1];

    double ref_x_prev = previous_path_x_[previous_path_x_.size() - 2];
    double ref_y_prev = previous_path_y_[previous_path_y_.size() - 2];

    ref_yaw = atan2(ref_y - ref_y_prev, ref_x-ref_x_prev);

    ptx.push_back(ref_x_prev);
    pty.push_back(ref_y_prev);

    ptx.push_back(ref_x);
    pty.push_back(ref_y);

  }

  // generate anchor points(from the current position of the car)
  vector<double> wp_1 = getXY(car_s_+30, getDforLane(curr_lane_), map_waypoints_s_, map_waypoints_x_, map_waypoints_y_);
  vector<double> wp_2 = getXY(car_s_+60, getDforLane(curr_lane_), map_waypoints_s_, map_waypoints_x_, map_waypoints_y_);
  vector<double> wp_3 = getXY(car_s_+90, getDforLane(curr_lane_), map_waypoints_s_, map_waypoints_x_, map_waypoints_y_);

  // push way points to the points
  ptx.push_back(wp_1[0]);
  pty.push_back(wp_1[1]);

  ptx.push_back(wp_2[0]);
  pty.push_back(wp_2[1]);

  ptx.push_back(wp_3[0]);
  pty.push_back(wp_3[1]);

  // transform all the way points to the ego car local co ordinate frame
  for (int i=0; i<ptx.size(); i++)
  {
    // sift reference frmae
    double shift_x = ptx[i] - ref_x;
    double shift_y = pty[i] - ref_y;

    ptx[i] = shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw);
    pty[i] = shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw);
  }

  // fit the spline among these waypoints
  tk::spline s;
  //std::sort(ptx.begin(), ptx.end());
  s.set_points(ptx, pty);

  // use the remaining points in the previous path
  // add these points to the new points, for smooth transition
  for (int i=0; i<previous_path_x_.size(); i++)
  {
    next_x_vals.push_back(previous_path_x_[i]);
    next_y_vals.push_back(previous_path_y_[i]);
  }

  // get the points from the spline and add to way points
  // set a horizon and get the points on segments for desired velocity
  double horizin_x = HORIZON_X_THRESHOLD;
  double horizin_y = s(horizin_x);
  double horizin_dist = distance(0, 0, horizin_x, horizin_y);

  double inc = 0;
  //  std::cout << "new points " << NO_OF_POINTS_PER_PATH-previous_path_x_.size() << "\n";
  // fill the rest of the points from the spline
  for (int i=0; i<NO_OF_POINTS_PER_PATH-previous_path_x_.size(); i++)
  {
    double N = horizin_dist/(dt_ * curr_vel_/MetersPerSec_To_MPH);
    double x_spline = inc + (horizin_x)/N;
    double y_spline = s(x_spline);

    inc = x_spline;

    // tansform the points back to map coordinates
    double x_ref = x_spline;
    double y_ref = y_spline;

    x_spline = x_ref*cos(ref_yaw)-y_ref*sin(ref_yaw);
    y_spline = x_ref*sin(ref_yaw)+y_ref*cos(ref_yaw);

    x_spline += ref_x;
    y_spline += ref_y;

    next_x_vals.push_back(x_spline);
    next_y_vals.push_back(y_spline);
  }

  // std::cout << "prev size: " << previous_path_x_.size() << " next size: " << next_x_vals.size() << "\n";

  //  std::cout << "final way points: \n";
  //  for (int i=0; i<next_x_vals.size(); i++)
  //  {
  //    std::cout << "x: " << next_x_vals[i] << " y: " << next_y_vals[i] << "\n";
  //  }
}
