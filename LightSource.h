// Tyrus Malmstrom
// Header file for the LightSource.cpp

#ifndef LIGHTSOURCE_H_INCLUDE
#define LIGHTSOURCE_H_INCLUDE

// directives:
#include <iostream>
#include <Eigen/Dense>
#include "Color.h"

// Namespace:
using namespace std;
using Eigen::Vector4d;

class LightSource {

 public:

  // Constructor:
  // Going to take a position (point), and an amount of Color :: red, green blue
 LightSource( const Vector4d&  _position, const Color&  _energy ):
  position( _position )
  ,energy  ( _energy   )
  { ++counter; };
  
  
  // pprint member function:
  void pprint(ostream& out = cout) const;


  // copy assignment operator: 1 of the BIG THREE
  const LightSource& operator= (const LightSource& rhs){
    if( this != &rhs ){ // Standard alias test...
      position = rhs.position;
      energy   = rhs.energy;
    }
    return *this;
  }
  
  
  Vector4d position;
  Color energy; 

 protected:
  static int counter; // to count the light sources in the scene.

  
};

//======================================================================

// output stream overloading:
ostream& operator<< (ostream& out, const LightSource& l);

#endif // LIGHTSOURCE_H_INCLUDE
