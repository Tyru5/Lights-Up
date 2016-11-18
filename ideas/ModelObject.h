// Author: Tyrus Malmstrom
// Date  : 11/1/2016
// Header file for pa4 ModelObject.cpp

#ifndef MODELOBJECT_H_INCLUDE
#define MODELOBJECT_H_INCLUDE


// directives:
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>

// namespace:
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ModelObject{

 public:
  // construtor:
  ModelObject(){};
  
  // accessors + mutators:
  int get_verticies()const;
  int get_faces()    const;
  double get_min_x() const;
  double get_max_x() const;
  double get_min_y() const;
  double get_max_y() const;
  double get_min_z() const;
  double get_max_z() const;

  VectorXd get_xes();
  VectorXd get_whys();
  VectorXd get_zeezs();
  
  MatrixXd get_main_vertex_list() const;

  // =======================
  void set_vertex_list(const MatrixXd& val);
  void set_faces_list(const MatrixXd& val);

  void set_verticies(const int& val);
  void set_faces(const int& val);
  // =======================
  
  // member functions:

  void print_vertex_list() const;

  double mean_vertex_x();
  double mean_vertex_y();
  double mean_vertex_z();

  void find_max_min_x();
  void find_max_min_y();
  void find_max_min_z();

  void print_faces_list() const;
  void extract_x_verts();
  void extract_y_verts();
  void extract_z_verts();

  double std_dev(const VectorXd& list);
  
  // class instance variables:
 private:
  int obj_verticies;
  int obj_faces;
  
  double total_x; /*Holds the total sum of all the x,y,z verticies*/
  double total_y;
  double total_z;

  double min_x;
  double max_x;

  double min_y;
  double max_y;

  double min_z;
  double max_z;

  VectorXd xes;
  VectorXd whys;
  VectorXd zeezs;

  MatrixXd vertex_list;
  MatrixXd faces_list;

};

#endif //MODELOBJECT_H_INCLUDE
