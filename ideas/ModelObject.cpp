// PA4 Assignment
// Author: Tyrus Malmstrom
// Date  : 11/1/2016
// Class : CS410
// Email : tyrus.alexander.malmstrom@gmail.com

// including directives:
#include <iostream>
#include <string>
#include <math.h> // for sqrt function
#include <vector>
#include <Eigen/Dense>
#include "ModelObject.h"
#include "Face.h"

// namespace
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int ModelObject:: get_verticies() const{
  return obj_verticies;
}

int ModelObject:: get_faces() const{
  return obj_faces;
}

void ModelObject:: set_verticies(const int& val){
  obj_verticies = val;
}

void ModelObject:: set_faces(const int& val){
  obj_faces = val;
}

void ModelObject::set_vertex_list(const MatrixXd& val){
  vertex_list = val;
}

void ModelObject::set_faces_list(const MatrixXd& val){
  faces_list = val;
}

double ModelObject:: mean_vertex_x(){
  total_x = 0.0;
  for (int i = 0; i < vertex_list.rows(); i++){
    total_x += vertex_list(i,0);
  }
  // cout << "Total_x = " << total_x << endl;
  double res = (total_x/obj_verticies);
  return res;
}

double ModelObject:: mean_vertex_y(){
  total_y = 0.0;
  for (int i = 0; i < vertex_list.rows(); i++){
    total_y += vertex_list(i,1);
  }
  // cout << "Total_y = " << total_y << endl;
  double res = (total_y/obj_verticies);
  return res;
}

double ModelObject:: mean_vertex_z(){
  total_z = 0.0;
  for (int i = 0; i < vertex_list.rows(); i++){
    total_z += vertex_list(i,2);
  }
  // cout << "Total_z = " << total_z << endl;
  double res = (total_z/obj_verticies);
  return res;
}

void ModelObject::find_max_min_x(){

  max_x = xes(0);
  min_x = xes(0);
  for (int i = 0; i < xes.size(); i++){
    if(xes(i) > max_x) max_x = xes(i);
    if(xes(i) < min_x) min_x = xes(i);
  }

  // cout << "This is the min x: " << min_x << endl;
  // cout << "This is the max x: " << max_x << endl;

}

double ModelObject:: get_min_x() const{
  return min_x;
}

double ModelObject:: get_max_x() const{
  return max_x;
}

void ModelObject::find_max_min_y(){

  max_y = whys(0);
  min_y = whys(0);
  for (int i = 0; i < whys.size(); i++){
    if(whys(i) > max_y) max_y = whys(i);
    if(whys(i) < min_y) min_y = whys(i);
  }

  // cout << "This is the min y: " << min_y << endl;
  // cout << "This is the max y: " << max_y << endl;

}

double ModelObject:: get_min_y() const{
  return min_y;
}

double ModelObject:: get_max_y() const{
  return max_y;
}

void ModelObject::find_max_min_z(){

  max_z = zeezs(0);
  min_z = zeezs(0);
  for (int i = 0; i < zeezs.size(); i++){
    if(zeezs(i) > max_z) max_z = zeezs(i);
    if(zeezs(i) < min_z) min_z = zeezs(i);
  }

  // cout << "This is the min z: " << min_z << endl;
  // cout << "This is the max z: " << max_z << endl;

}

double ModelObject:: get_min_z() const{
  return min_z;
}

double ModelObject:: get_max_z() const{
  return max_z;
}

void ModelObject:: print_vertex_list() const{
  cout << vertex_list << endl;
}

void ModelObject::print_faces_list() const{
  cout << faces_list << endl;
}

void ModelObject::extract_x_verts(){

  xes.resize( obj_verticies );
  for (int i = 0; i < vertex_list.rows(); i++){
    xes(i) = vertex_list(i,0);
  }

}

void ModelObject::extract_y_verts(){

  whys.resize( obj_verticies );
  for (int i = 0; i < vertex_list.rows(); i++){
    whys(i) = vertex_list(i,1);
  }

}


void ModelObject::extract_z_verts(){

  zeezs.resize(obj_verticies);
  for (int i = 0; i < vertex_list.rows() ; i++){
    zeezs(i) = vertex_list(i,2);
  }

}

VectorXd ModelObject:: get_xes(){
  return xes;
}

VectorXd ModelObject:: get_whys(){
  return whys;
}

VectorXd  ModelObject:: get_zeezs(){
  return zeezs;
}


double ModelObject::std_dev(const VectorXd& list){
  // will hold the total value of the list.
  double total_sum = 0;
  for(int i = 0; i < obj_verticies; i++) {
    total_sum += list(i);
  }

  // calculating the mean:
  double mean = total_sum / obj_verticies;

  double element_mean_squared = 0.0;
  for(int i = 0; i < obj_verticies; i++){
    element_mean_squared += pow( (list(i) - mean), 2 );
  }

  double res = (element_mean_squared / obj_verticies);

  return sqrt(res);
}

MatrixXd ModelObject:: get_main_vertex_list() const{
  return vertex_list;
}

