#include <traj_planner.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"

trajPlanner::trajPlanner()
{
 std::cout << "trajectory planner constrcuted \n";
}

trajPlanner::~trajPlanner()
{

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

void trajPlanner::generateTrajctory(double car_x, double car_y, double car_yaw, std::vector<double>& next_x_vals, std::vector<double>& next_y_vals)
{

    double dist_inc = 0.5;
    for(int i = 0; i < 50; i++)
    {
          next_x_vals.push_back(car_x + (dist_inc*i)*std::cos(car_yaw*M_PI/180));
          next_y_vals.push_back(car_y + (dist_inc*i)*std::sin(car_yaw*M_PI/180));
    }

    std::cout << "final way points: \n";
    for (int i=0; i<next_x_vals.size(); i++)
    {
        std::cout << "x: " << next_x_vals[i] << " y: " << next_y_vals[i] << "\n";
    }
}

