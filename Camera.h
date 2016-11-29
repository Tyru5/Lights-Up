// Tyrus Malmstrom
// Header file for the Camera.cpp

#ifndef CAMERA_H_INCLUDE
#define CAMERA_H_INCLUDE

// directives:
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <tuple> 
#include <Eigen/Dense>
#include <Eigen/Geometry> // for cross product of vectors.
#include "ModelObject.h"
#include "Ray.h"
#include "Face.h"
#include "Color.h"
#include "LightSource.h"
#include "Sphere.h"


using Eigen::Vector3d;
using Eigen::Vector3i;
using Eigen::RowVector3i;

class Camera{

 public:

  // Defining some list types to hold multiple:
  typedef vector<ModelObject> ModelObjects;
  typedef vector<LightSource> LightSources;
  typedef vector<Sphere> SphereList;

  // constructor:
  Camera(){};

  // member functions:
  void parseScene( const string& scene_file );
  void buildRM();
  void calculateRays();

  void print_ts(const vector< vector<double>>& vect);
  void find_tmin_tmax( vector< vector<double>>& tvals);
  
  // Where the magic happens:  
  RowVector3i mapColour(const Color& bc);
  void writeSpheresAndModels( const string& out_file );
  void writeSpheres(const string& out_file);
  void writeModels( const string& out_file );
  void printPixs() const;
  
  void rayTriangleIntersection();
  void computeDist( const Face& current_face );

  void print_ptof();
  
  // Methods for adding objects:~~~~~~~~~~~~~~~~~~~~~
  void addLightSources( const LightSource& light ){
    lightSource_list.push_back( light );
  }
  
  void addModels( const ModelObject& model ){
    modelObject_list.push_back( model );
  }
  
  void addSphere( const Sphere& sphere ){
    spheres.push_back( sphere );
  }
  // done for adding objects~~~~~~~~~~~~~~~~~~~~~~~~

  const int numberOfSpheres() const{
    return static_cast<int>(spheres.size());
  }

  const int numberOfModels() const{
    return static_cast<int>(modelObject_list.size());
  }
  
  
  // class instance variables:
 protected:
  // location of the focal point
  string eye_header;
  // the look at point
  string look_header;
  // up vector
  string supv;
  // distacne from ip:
  string dist_header;
  // bounds
  string bounds_header;
  // res
  string res_header;

  // Camera specs:
  Vector3d EYE;
  Vector3d LOOKAP;
  Vector3d UPV;

  // RM basis vectors:
  Vector3d WV;
  Vector3d UV;
  Vector3d VV;

  double dist;

  double bottom;
  double left;
  double top;
  double right;

  double width; // had to change this to double to have correct division
  double height;

  // array of rays hey...
  vector< vector< Ray > > Rays;

  // 2d array to hold all t's:
  vector< vector< double > > ts;
  vector< vector< Color > > ptof;

  vector< vector< RowVector3i> > sphere_pixs;
  vector< vector< RowVector3i> > model_pixs;
  vector< vector< RowVector3i> > sphere_model_pixs;

  double tmin = numeric_limits<double>::max(); // max double
  double tmax = numeric_limits<double>::min(); // min double

  // ambient illumination in the scene:
  Color ambient_color;

  // Actual Lists:
  ModelObjects modelObject_list;
  LightSources lightSource_list;
  SphereList spheres;
  
};

#endif // CAMERA_H_INCLUDE
