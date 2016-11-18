// Tyrus Malmstrom
// 11/9/2016
// Color.cpp class for handling and manipulating color


// directives:
#include <iostream>
#include "Color.h"

// namespace:
using namespace std;


// Macros:
#define DEBUG false


void Color::pprint(ostream& out) const{
  out << "(" << red << ", " << green << ", " << blue << ")" << endl; 
}

ostream& operator<< (ostream& out, const Color& c){
  c.pprint( out );
  return out;
}
