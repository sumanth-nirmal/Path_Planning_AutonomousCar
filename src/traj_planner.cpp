#include <traj_planner.h>

trajPlanner::trajPlanner()
{

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

// source: http://mplab.ucsd.edu/tutorials/minimumJerk.pdf
// start and end is a vector with position, velocity and acceleration
void trajPlanner::minJerkTrajParam(double start[3], double end[3], double t, double& traj[6])
{
    // first 3 co-effiecients from the initial conditions
    double A0 = start[0];
    double A1 = start[1];
    double A2 = start[2]/2;

    // last 3 co - efficients from terminal conditions
    // minimum jerk in matrix form
    Eigen::Matrix3Xf T;

    T << std::pow(t, 3),      std::pow(t, 4),      std::pow(t, 5),
         3 * std::pow(t, 3),  4 * std::pow(t, 4),  5 * std::pow(t, 5),
         6 * t,               12 * std::pow(t, 2), 20 * std::pow(t, 3);

    Eigen::MatrixXf X(3,1);

    X << end[0] - A0 - A1*t - A2*t*t,
         end[1] - A1 - 2*A2*t,
         end[2] - 2*A2;

   Eigen::MatrixXf A = T.inverse() * X;

   // parameters for minimum jerk trajectory
   traj[0] = A0;   traj[1] = A1;     traj[2] = A2;
   traj[3] = X[0]; traj[4] = X[1];   traj[5] = X[3];

}

void trajPlanner::generateTrajctory(std::vector<double>& next_x_vals, std::vector<double>& next_y_vals)
{

}

