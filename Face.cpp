// Tyrus Malmstrom
// 11/1/2016
// Face.cpp class 


// directives:
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Face.h"

// namespace:
using namespace std;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;

// Macros:
#define DEBUG false

void Face::map(const MatrixXd& mat){

  mvil.resize(3,3);

  /*
  cout << "this is A \n" << A << endl;
  cout << "this is B \n" << B << endl;
  cout << "this is C \n" << C << endl;

  cout << "the A row in mat \n" << mat.row(A).transpose() << endl;
  cout << "the B row in mat \n" << mat.row(B).transpose() << endl;
  cout << "the C row in mat \n" << mat.row(C).transpose() << endl;
  */
  
  mvil.col(0) = mat.row(A).transpose();
  mvil.col(1) = mat.row(B).transpose();
  mvil.col(2) = mat.row(C).transpose();

  // cout << mvil << endl;

}

void Face::pprint(ostream& out) const{
  out << "FACE: " << endl;
  // print off vertex list:
  out << "mapped vertex list: " << endl;
  out << mvil << endl;
}


ostream& operator<< (ostream& out, const Face& f){
  f.pprint( out );
  return out;
}

Vector3d Face::getA() const{
  return mvil.col(0);
}

Vector3d Face::getB() const{
  return mvil.col(1);
}

Vector3d Face::getC() const{
  return mvil.col(2);
}

void Face::addFace(const Face& f){
  Faces.push_back(f);
}

Face Face::getFace(const int& index) const{
  return Faces[index];
}