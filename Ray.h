// Tyrus Malmstrom
// Header file for the Ray.cpp

#ifndef RAY_H_INCLUDE
#define RAY_H_INCLUDE

// directives:
#include <iostream>
#include <limits>
#include <vector>
#include <Eigen/Dense>


using std::cout;
using std::ostream;
using std::endl;
using std::vector;
using Eigen::Vector3d;

class Ray{

 public:
  // Defualt Constructor:
 Ray(): origin( Vector3d() ), direction( Vector3d() ) {};
 Ray(const Vector3d& origin_, const Vector3d& direction_):
  origin( origin_ ), direction( direction_ ) {};

  // pprint member function:
  void pprint(ostream& out = cout) const;
    
 // copy assignment operator: 1 of the BIG THREE
  const Ray& operator= (const Ray& rhs){
    if( this != &rhs ){ // Standard alias test...
      origin = rhs.origin;
      direction = rhs.direction;
    }
    return *this;
  }

  Vector3d origin;
  Vector3d direction;

};

ostream& operator<< (ostream& out, const Ray& r);

#endif // RAY_H_INCLUDE
