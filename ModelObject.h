// Author: Tyrus Malmstrom
// Date  : 11/1/2016
// Header file for pa4 ModelObject.cpp

#ifndef MODELOBJECT_H_INCLUDE
#define MODELOBJECT_H_INCLUDE

// directives:
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "Face.h"

/* Using tinyobjloader from syoyo*/
#include "tiny_obj_loader.h"


using std::string;
using std::vector;
using Eigen::MatrixXd;

class ModelObject{

 public:
  // construtor:
 ModelObject( const string& _obj_file, const double _tx, const double& _ty, const double& _tz, const double& _wx, const double& _wy, const double& _wz, const double& _theta ):
  obj_file( _obj_file )
  ,tx( _tx )
  ,ty( _ty )
  ,tz( _tz )
  ,wx( _wx )
  ,wy( _wy )
  ,wz( _wz )
  ,theta( _theta )
  {};

  
  // pprint member function:
  void pprint(ostream& out = cout) const;

  // Member functions:
  void parseObj();
  void PrintInfo()const;
  void getFaces();
  void rayTriangleIntersection(const int& width, const int& height);
  void computeDist();
  void print_ts(const vector<vector<double>>& vect);

  
 protected:
  /* For parsing obj file(s)*/
  string obj_file;
  tinyobj::attrib_t attrib;
  vector<tinyobj::shape_t> shapes;
  vector<tinyobj::material_t> materials;
  string err;


  /* For translatoin / transformation and axis angle rotation */
  double tx;
  double ty;
  double tz;

  double wx;
  double wy;
  double wz;
  double theta;

  MatrixXd vertices;
  vector<Face> F;

  // 2d array to hold all t's:
  vector<vector< double > > ts; 

};


// output stream overloading:
ostream& operator<< (ostream& out, const ModelObject& m);


#endif //MODELOBJECT_H_INCLUDE
