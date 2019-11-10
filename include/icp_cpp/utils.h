#include <cmath>
#include <eigen3/Eigen/Dense>


double degrees(double radians)
{
  long double pi = 3.141592653589793238463;
  return radians * 180. / pi;
}

double radians(double degrees)
{
  long double pi = 3.141592653589793238463;
  return degrees * pi/180.;
}


Eigen::Vector3d matrix_to_rpy(Eigen::MatrixXd rot_matrix, bool degrees_flag=false)
{
  double roll = std::atan(rot_matrix(2,1)/rot_matrix(2,2));
  double pitch = std::atan(-rot_matrix(2,0)/sqrt(rot_matrix(2,1)*rot_matrix(2,1) + rot_matrix(2,2)*rot_matrix(2,2)) );
  double yaw = std::atan(rot_matrix(1,0)/rot_matrix(0,0));

  if (degrees)
  {
    roll = degrees(roll);
    pitch = degrees(pitch);
    yaw = degrees(yaw);

  }
  Eigen::Vector3d out(roll, pitch, yaw);
  return out;

}
