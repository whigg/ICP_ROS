#include <eigen3/Eigen/Dense>
#include <vector>
#include <iostream>
#include "icp_cpp/icp.h"
#include "icp_cpp/utils.h"
#include "nanoflann/nanoflann.hpp"
#include <numeric>
#include <cstdlib>
#include <typeinfo>
#include <cmath>




ICP::ICP(int m_it)

{
  max_iters = m_it;
  indices = {};
  dists = {};

}


void print_mat_dims(const Eigen::MatrixXd &m1, const Eigen::MatrixXd &m2, const std::string&& m1_name, const std::string&& m2_name)
{
  std::cout<<m1_name<<" is ("<<m1.rows()<<", "<<m1.cols()<<") | "<<m2_name<<" is ("<<m2.rows()<<", "<<m2.cols()<<")"<<std::endl;

}
Eigen::MatrixXd ICP::best_fit_transform(const Eigen::MatrixXd &pointcloud_A, const Eigen::MatrixXd &pointcloud_B)
{
  Eigen::MatrixXd hg_T = Eigen::MatrixXd::Identity(4,4);

  std::cout<<"A Trace: "<<pointcloud_A.trace()<<std::endl;
  std::cout<<"B Trace: "<<pointcloud_B.trace()<<std::endl;
  Eigen::MatrixXd shifted_A = pointcloud_A;
  Eigen::MatrixXd shifted_B = pointcloud_B;
  Eigen::Vector3d cent_A(0.,0.,0.), cent_B(0.,0.,0.);
  std::cout<<"S A Trace: "<<shifted_A.trace()<<std::endl;
  std::cout<<"S B Trace: "<<shifted_B.trace()<<std::endl;

  int a_rows = pointcloud_A.rows(), a_cols = pointcloud_A.cols(), b_rows=pointcloud_B.rows(), b_cols = pointcloud_B.cols();

  std::cout<<"Dims ("<<a_rows<<", "<<a_cols<<") and ("<<b_rows<<", "<<b_cols<<") "<<std::endl;
  for (int i = 0; i < a_rows; i++)
  {
    cent_A += pointcloud_A.block<1,3>(i,0).transpose();
    cent_B += pointcloud_B.block<1,3>(i,0).transpose();
  }

  cent_A/=a_rows;
  cent_B/=a_rows;

  std::cout<<"Centroid A: "<< cent_A.transpose()<<std::endl;
  std::cout<<"Centroid B: "<< cent_B.transpose()<<std::endl;

  for(int i = 0; i < a_rows; i++)
  {
    shifted_A.block<1,3>(i,0) = pointcloud_A.block<1,3>(i,0) - cent_A.transpose();
    shifted_B.block<1,3>(i,0) = pointcloud_B.block<1,3>(i,0) - cent_B.transpose();
  }


  Eigen::MatrixXd ata = shifted_A.transpose()*shifted_A;

  Eigen::MatrixXd btb = shifted_B.transpose()*shifted_B;

  Eigen::MatrixXd H = shifted_A.transpose()*shifted_B;


  std::cout<<"Resulting H: "<<std::endl<<H<<std::endl;

  Eigen::MatrixXd U;
  Eigen::VectorXd S;
  Eigen::MatrixXd V;
  Eigen::MatrixXd Vt;
  Eigen::MatrixXd rot_matrix;
  Eigen::Vector3d trans_vector;

  Eigen::JacobiSVD<Eigen::MatrixXd> SVD(H, Eigen::ComputeFullU | Eigen::ComputeFullV);

  U = SVD.matrixU();
  S = SVD.singularValues();
  V = SVD.matrixV();
  Vt = V.transpose();

  std::cout<<"U: \n"<<U<<" \n S: \n" <<S<<"\nV: \n"<<V<<std::endl;

  rot_matrix = V*U.transpose();

  Eigen::Vector3d rpy = matrix_to_rpy(rot_matrix);

  rpys.push_back(rpy);

  std::cout<<"Rotation Vector: "<<rpy.transpose()<<std::endl<<std::endl;
  std::cout<<"Rotation: \n \n"<<rot_matrix<<std::endl;
  if (rot_matrix.determinant() < 0)
  {

    Vt.block<1,3>(2,0) *=-1;

    rot_matrix = Vt*U.transpose();


  }

  trans_vector = cent_B - rot_matrix * cent_A;
  hg_T.block<3,3>(0,0) = rot_matrix;
  hg_T.block<3,1>(0,3) = trans_vector;

  return hg_T;
}

float ICP::euc_dist(const Eigen::Vector3d &pt_a, const Eigen::Vector3d &pt_b)
{
  return sqrt(pow(pt_a[0]-pt_b[0], 2) + pow(pt_a[1]-pt_b[1], 2) + pow(pt_a[2]-pt_b[2], 2) );
}

void ICP::calc_closest_neighbors_kdtree(const Eigen::MatrixXd &src, const nn_kd_tree &dst_tree)
{
  const size_t num_results = 1;
  Eigen::Vector3d src_pt;

  for (int i = 0; i < src.rows(); i++)
  {
    std::vector<size_t> temp_index(1);
    std::vector<double> temp_dist(1);


    src_pt = src.block<1,3>(i,0).transpose();
    nanoflann::KNNResultSet<double> results(num_results);
    results.init(&temp_index[0], &temp_dist[0]);

    std::vector<double> query_pt{src_pt(0), src_pt(1), src_pt(2)};
    dst_tree.index->findNeighbors(results, &query_pt[0], nanoflann::SearchParams(10));

    indices.push_back(int(temp_index[0]));
    dists.push_back(temp_dist[0]);
  }

}



void ICP::run_scan_matcher(const Eigen::MatrixXd &pointcloud_A, const Eigen::MatrixXd &pointcloud_B, double tolerance)
{

  rpys.clear();

  int a_rows = pointcloud_A.rows();
  int b_rows = pointcloud_B.rows();
  Eigen::MatrixXd hg_src = Eigen::MatrixXd::Ones(3+1, a_rows);
  src_3d = Eigen::MatrixXd::Ones(3, a_rows);
  Eigen::MatrixXd hg_dst = Eigen::MatrixXd::Ones(3+1, b_rows);

  Eigen::MatrixXd closest_pts_in_dst = Eigen::MatrixXd::Ones(3,a_rows);

  const size_t dim = 3;

  std::cout<<"Building KD Tree"<<std::endl;
  nn_kd_tree dst_index(dim, std::cref(pointcloud_B), 10);
  dst_index.index->buildIndex();

  std::cout<<"Built KD Tree"<<std::endl;



  print_mat_dims(hg_src, pointcloud_A, "Hg_Src","Input PC_A");


  print_mat_dims(src_3d, src_3d, "src_3d", "src_3d again");
  for (int i = 0; i < a_rows; i++)
  {
    //std::cout<<"I is "<<i<<std::endl;

    hg_src.block<3,1>(0,i) = pointcloud_A.block<1,3>(i,0).transpose();
    src_3d.block<3,1>(0,i) = pointcloud_A.block<1,3>(i,0).transpose();


  }

  print_mat_dims(hg_dst, pointcloud_B, "Hg_Dst", "Input PC_B");
  for (int j = 0; j < b_rows; j++)
  {
    hg_dst.block<3,1>(0,j) = pointcloud_B.block<1,3>(j,0).transpose();
  }



  std::cout<<"Calculating closest neighbors..."<<std::endl;
  double prev_error = 0.0;
  double mean_error = 0.0;

  double diff;

  for (iters = 1; iters <= max_iters; iters++)
  {
    indices.clear();
    dists.clear();
    calc_closest_neighbors_kdtree(src_3d.transpose(), dst_index);
    std::cout<<"Finished calculating closest. "<<std::endl;

    print_mat_dims(hg_dst, closest_pts_in_dst, "Hg_Dst", "Clos_Index");
    for (int j = 0; j < a_rows; j++)
    {
        closest_pts_in_dst.block<3,1>(0,j) = hg_dst.block<3,1>(0,indices[j]);
    }




    std::cout<<"Doing best fit transform..."<<std::endl;
    transform_matr = best_fit_transform(src_3d.transpose(), closest_pts_in_dst.transpose());


    hg_src = transform_matr*hg_src;

    for (int j = 0; j < a_rows; j++)
      src_3d.block<3,1>(0,j) = hg_src.block<3,1>(0,j);

    mean_error = std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
    diff = mean_error-prev_error;


    diff = (diff < 0 ) ? diff*-1 : diff;

    std::cout<<"Prev Error: "<<prev_error<<", "<<typeid(prev_error).name()<<", Mean Error: "<<mean_error<<", "<<typeid(mean_error).name()<<", Diff: "<<diff<<", "<<typeid(diff).name()<<", Tolerance: "<<tolerance<<std::endl;
    if (diff < tolerance)
      break;

    prev_error = mean_error;

    iters++;

  }
  std::cout<<"Final Prev Error:"<<prev_error<<std::endl;
  std::cout<<"Final Mean Error: "<<mean_error<<std::endl;
  std::cout<<"Final Diff Error: "<<diff<<std::endl;

  for (auto rpy : rpys)
    std::cout<<rpy;
  std::cout<<std::endl;
  transform_matr = best_fit_transform(pointcloud_A, src_3d.transpose());
}
