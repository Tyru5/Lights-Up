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
#include "ModelObject.h"
#include "Ray.h"
#include "Face.h"
#include "Color.h"
#include "LightSource.h"
#include "Sphere.h"

using Eigen::Vector3d;
using Eigen::Vector3i;

class Camera{

 public:

  // Defining some list types to hold multiple:
  typedef vector<ModelObject> ModelObjects;
  typedef vector<LightSource> LightSources;
  typedef vector<Sphere> SphereList;

  // constructor:
  Camera(){};

  // member functions:
  void parseScene( const std::string& scene_file );

  void buildRM();
  void calculateRays();
  
  void print_ts(const std::vector<std::vector<double>>& vect);
  void find_tmin_tmax(std::vector<std::vector<double>>& tvals);
  
  // Where the magic happens:
  void computeDist(const Face& current_face);
  
  Vector3i getColour(const double& tval);
  Vector3i mapColour(const std::tuple<bool, Color>& bc);
  
  void writeImage(const std::string& out_file);
  
  void rayTriangleIntersection(const ModelObject& obj, const Face& face);

  void raySphereIntersection();
  

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

  
  // class instance variables:
 protected:
  // location of the focal point
  std::string eye_header;
  // the look at point
  std::string look_header;
  // up vector
  std::string supv;
  // distacne from ip:
  std::string dist_header;
  // bounds
  std::string bounds_header;
  // res
  std::string res_header;

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
  std::vector< std::vector< Ray > > Rays;

  // 2d array to hold all t's:
  std::vector< std::vector< double > > ts; 


  // TESTING:
  std::vector< tuple<bool, Color> > sphere_colors;
  
  double tmin = std::numeric_limits<double>::max(); // max double
  double tmax = std::numeric_limits<double>::min(); // min double

  // ambient illumination in the scene:
  Color ambient_color;

  // Actual Lists:
  ModelObjects modelObject_list;
  LightSources lightSource_list;
  SphereList spheres;
  
};

#endif // CAMERA_H_INCLUDE
