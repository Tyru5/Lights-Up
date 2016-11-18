// Tyrus Malmstrom
// Header file for the Face.cpp

#ifndef FACE_H_INCLUDE
#define FACE_H_INCLUDE

// directives:
#include <iostream>
#include <string>
#include <Eigen/Dense>

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;

// for overloading operator<<
using std::cout;
using std::ostream;

class Face{

 public:
  // Default Constructor:
 Face(): A(0),B(0),C(0){};
 Face(const double& A_, const double& B_, const double& C_): A(A_),B(B_),C(C_){};


  // Member functions:
  void map(const MatrixXd& mat);
  void pprint(ostream& out = cout) const;
  void addFace(const Face& f);
  Face getFace(const int& index) const;
  
  Vector3d getA() const;
  Vector3d getB() const;
  Vector3d getC() const;
  
  // copy assignment operator: 1 of the BIG THREE
  // This doesn't really make sense yet...
  const Face& operator= (const Face& rhs){
    if( this != &rhs ){ // Standard alias test...
      A = rhs.A;
      B = rhs.B;
      C = rhs.C;
      mvil = rhs.mvil;
    }
    return *this;
  }

 public:
  double A;
  double B;
  double C;

  Matrix3d mvil;

  std::vector< Face > Faces;

};

ostream& operator<< (ostream& out, const Face& f);

#endif // FACE_H_INCLUDE
