// Tyrus Malmstrom
// Header file for the Color.cpp

#ifndef COLOR_H_INCLUDE
#define COLOR_H_INCLUDE

// directives:
#include <iostream>

// Namespace:
using namespace std;


class Color {

 public:

  // class instance variables:
  // colors values represented as a double :: colors for this assignment should be in the range (inclusive) of 0-1:
  double red;
  double green;
  double blue;

  
  // Constructor(s):
 Color():
  red( 0.0 )
    ,green( 0.0 )
    ,blue( 0.0 )
    {};

 Color(const double& _red, const double& _green, const double& _blue):
  red(    _red )
    ,green( _green )
    ,blue(  _blue )
  {};


  // pprint member function:
  void pprint(ostream& out = cout) const;

  
  // operator overloading:
  Color& operator += (const Color& other_color){

    red += other_color.red;
    green += other_color.green;
    blue +=  other_color.blue;
    return *this;
    
  }

  Color& operator += (const double& factor){

    red   += factor;
    green += factor;
    blue  += factor;
    return *this;
    
  }


  
  Color& operator *= (const Color& other_color){
    
    red *= other_color.red;
    green *= other_color.green;
    blue *= other_color.blue;
    return *this;
    
  }


  Color& operator *= (const double& factor){
    
    red *= factor;
    green *= factor;
    blue *= factor;
    return *this;
    
  }

  

  Color& operator /= (const double& num){

    red /= num;
    green /= num;
    blue /= num;
    return *this;
    
  }
  
  // copy assignment operator: 1 of the BIG THREE
  const Color& operator= (const Color& rhs){
    if( this != &rhs ){ // Standard alias test...
      
      red  = rhs.red;
      blue = rhs.blue;
      green= rhs.green;
      
    }
    return *this;
  }


  // member function to check validity of color values:
  void validate_colors() const{
    if( red < 0.0 || green < 0.0 || blue < 0.0 ){
      cerr << "This is not a valid color!" << endl;
    }
  }


  
};

//======================================================================

inline Color operator * (const Color& f_color, const Color& s_color){

  return Color(
	       f_color.red   *   s_color.red,
	       f_color.green * s_color.green,
	       f_color.blue  *  s_color.blue 
	       );
    
}

inline Color operator + (const Color& f_color, const Color& s_color){

  return Color(
	       f_color.red   +   s_color.red,
	       f_color.green + s_color.green,
	       f_color.blue  +  s_color.blue 
	       );
    
}

inline Color operator * ( const Color& f_color, const double& scalar){

  return Color(
	       f_color.red   *   scalar,
	       f_color.green * scalar,
	       f_color.blue  *  scalar
	       );
    
}


// output stream overloading:
ostream& operator<< (ostream& out, const Color& c);

#endif // COLOR_H_INCLUDE
