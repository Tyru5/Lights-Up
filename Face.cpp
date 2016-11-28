// Tyrus Malmstrom
// 11/1/2016
// Face.cpp class 


// directives:
#include "Face.h"

// namespace:
using namespace std;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;

// Macros:
#define DEBUG false

void Face::map(const MatrixXd& mat, const Matrix3d& fm, const Vector3d& sn){

  mvil.resize(3,3);
  material.resize(3,3);

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

  // mapping face material:
  material = fm;
  // mapping surface normal:
  surface_normal = sn;
}

void Face::pprint(ostream& out) const{
  out << "FACE: " << endl;
  // print off vertex list:
  out << "mapped vertex list: " << endl;
  out << mvil << endl;
  out << "With this associated material: \n" << material << endl;
  out << "With surface normal: \n" << surface_normal << endl;
  out << "With ptof of = \n" << ptof << endl;
}


ostream& operator<< (ostream& out, const Face& f){
  f.pprint( out );
  return out;
}

const Vector3d Face::getA() const{
  return mvil.col(0);
}

const Vector3d Face::getB() const{
  return mvil.col(1);
}

const Vector3d Face::getC() const{
  return mvil.col(2);
}
