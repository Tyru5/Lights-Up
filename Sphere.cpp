// Tyrus Malmstrom
// 11/9/2016
// Sphere.cpp class for handling Spheres


// directives:
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <Eigen/Dense>
#include "Sphere.h"

// namespace:
using namespace std;
using Eigen::Vector3d;

// Macros:
#define DEBUG true


int Sphere::counter = 0;

void Sphere::pprint(ostream& out) const{
  out << "~Sphere " << counter << "~\ncenter : " << center.transpose() << "\nradius = " << radius << "\ncolor = " << mat_props;
}

ostream& operator<< (ostream& out, const Sphere& s){
  s.pprint( out );
  return out;
}


bool Sphere::raySphereIntersection( const Ray& ray ){

  /*
    Using the Faster Method algorithm presented in class
    Returns the closest intersection point by defualt.
  */

  Vector3d ray_origin = ray.origin;
  Vector3d ray_direction = ray.direction;
  Vector3d Tv = center - ray_origin;
  // if(DEBUG) cout << "base to center = \n" << Tv << endl;
  double v    = Tv.dot( ray_direction );
  double csq  = Tv.dot(Tv);
  double disc = ( pow( radius, 2.0 ) - (csq - pow(v, 2.0) ) );
  // if(DEBUG) cout << "v,csq,disc = " << v << " " << csq << " " << disc << endl;
  if( disc < 0.0 ){
    return false;
  }
  else{
    d = sqrt(disc);
    t = v - d;
    // if(DEBUG) cout << "tval = " << t << endl;
    ptos = ray_origin + t * ray_direction;
    // cout << "intersection point =" << ptos.transpose() << endl;
    
  }

  return true;
  
}

tuple<bool, Color> Sphere::getRaySphereRGB( const Ray& ray, const Color& ambl, const vector<LightSource>& lights ){

  /*
    Given a certain ray-sphere intersection, compute the RGB off the surface:
  */

  bool res = raySphereIntersection( ray );
  double alpha = 16.0;
  Color color; // to start off, a blank color;
  if( res ){

    Vector3d snrm = ptos - center; snrm = snrm/snrm.norm(); // JUST UPDATED IT!
    // if(DEBUG) cout << "the snrm on sphere is = " << snrm.transpose() << " with ptos = " << ptos.transpose() << endl;
    // Initial condition of the ambient lighting of the scene:
    color = ambl * mat_props;
    // cout << color;

    // cout << lights.size() << endl;
    for( int z = 0; z < static_cast<int>( lights.size() ); z++){
    
      Vector3d lp( lights[z].position(0), lights[z].position(1), lights[z].position(2) );
      // if(DEBUG) cout << "light position = " << lp.transpose() << endl;
    
      Vector3d toL = lp - ptos; toL = toL/toL.norm(); // unit length
      // cout << "toL = " << toL.transpose() << " with associated ptos = " << ptos.transpose() << endl;
    
      if( snrm.dot( toL ) > 0.0 ){ // meaning there is actually an angle

	color += mat_props * lights[z].energy * snrm.dot( toL );
	// cout << "color2 = " << color;
	Vector3d toC  = ray.origin - ptos; toC = toC / toC.norm();
	// cout << "toC = " << toC.transpose() << " with associated ptos = " << ptos.transpose() << endl;
	
	Vector3d spR  = (2 * snrm.dot( toL ) * snrm) - toL;
	// cout << "spR = " << spR.transpose() << " with ptos of = " << ptos.transpose() << endl;

	// cout << toC.dot( spR ) << " ptos associated = " << ptos.transpose() << endl;; //<-- why not 16?
	
	color += mat_props * lights[z].energy *  pow( toC.dot( spR ), alpha );
	// cout << "color3 = " << color << "with ptos of = " << ptos.transpose() << endl;

      }

    }
    
    // cout << "about to return the color." << endl;
    return make_tuple(true, color);
    
  }else{

    return make_tuple(false, Color() );
    
  }

}
