// Tyrus Malmstrom
// 11/9/2016
// LightSource.cpp class for handling the light sources in the scene


// directives:
#include "LightSource.h"
#include <iostream>

// namespace:
using namespace std;


// Macros:
#define DEBUG false

int LightSource::counter = 0;

void LightSource::pprint(ostream& out) const{
  out << "~Light Source " << counter << "~" << endl;
  out << "Located at: \n" << position << endl;
  out << "\"Energy\" emmisison of red, green, blue: " << energy;
}

ostream& operator<< (ostream& out, const LightSource& l){
  l.pprint( out );
  return out;
}
