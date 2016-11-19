// Tyrus Malmstrom
// Header file for the Sphere.cpp

#ifndef SPHERE_H_INCLUDE
#define SPHERE_H_INCLUDE

// directives:
#include <iostream>
#include <string>
#include <tuple> // std::tuple, std::get, std::tie, std::ignore
#include <Eigen/Dense>
#include "Color.h"
#include "Ray.h"
#include "LightSource.h"

// namespace:
using Eigen::Vector3d;


using std::vector;

class Sphere {

 public:
  
  // Constructor(s):
 Sphere():
  center( Vector3d() )
 ,radius( 0.0 )
 ,mat_props( 0.0,0.0,0.0 )
 {counter++;};

 Sphere( const Vector3d& _center, const double& _radius, const Color& _mat_props ):
  center( _center )
 ,radius( _radius )
 ,mat_props( _mat_props )
 {counter++;};
  
  // Member functions:
  bool raySphereIntersection( const Ray& ray );
  tuple<bool, Color> getRaySphereRGB( const Ray& ray, const Color& ambl, const vector<LightSource>& lights );

  // pprint member function:
  void pprint(ostream& out = cout) const;

  
  const Vector3d getCenter() const{
    return center.transpose();
  }

  const double getRadius() const{
    return radius;
  }

  // copy assignment operator: 1 of the BIG THREE
  const Sphere& operator= (const Sphere& rhs){
    if( this != &rhs ){ // Standard alias test...
      center = rhs.center;
      radius = rhs.radius;
    }
    return *this;
  }


  
  // class instance variables:
 protected:
  Vector3d center;
  double radius;
  Color mat_props;

  double d;
  double t;
  Vector3d ptos;

  
  static int counter; // to count the spheres in the scene.
  
};

ostream& operator<< (ostream& out, const Sphere& s);

#endif // SPHERE_H_INCLUDE
