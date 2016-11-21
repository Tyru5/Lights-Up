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
#include "Face.h"

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

  const int getFaces() const{
    return number_of_faces;
  }

  const int getVerts() const{
    return number_of_vertices;
  }
  
  const int getVNorms() const{
    return number_of_vertices_norm;
  }
  
 protected:
  string obj_file;

  /* For translatoin / transformation and axis angle rotation */
  double tx;
  double ty;
  double tz;

  double wx;
  double wy;
  double wz;
  double theta;

  // v's
  int number_of_vertices;
  MatrixXd vertices;

  // vn's
  int number_of_vertices_norm;
  MatrixXd nertices;

  // faces -- not really sure yet:
  int number_of_faces;
  vector<Face> F;

  
};


// output stream overloading:
ostream& operator<< (ostream& out, const ModelObject& m);


#endif //MODELOBJECT_H_INCLUDE
