// Tyrus Malmstrom
// Header file for the Scene.cpp

#ifndef SCENE_H_INCLUDE
#define SCENE_H_INCLUDE

// directives:
#include <iostream>
#include <string>
#include <vector>
#include "Camera.h"
#include "PlyImage.h"
#include "ModelObject.h"
#include "Face.h"
#include "Color.h"
#include "LightSource.h"
#include "Sphere.h"

class Scene {

 public:

  // Defining some list types to hold multiple:
  typedef vector<ModelObject> ModelObjects;
  typedef vector<LightSource> LightSources;
  typedef vector<Sphere> SphereList;

  
  // Constructor:
  // Default background color is black
 Scene( const Color& _background_color = Color() ):
  background_color( _background_color )
  {};

  // Member functions:
  void parseScene( const string& scene_file, Camera& camera );

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

  
  const Color& get_ambient() const {
    return ambient_color;
  }

  // getters for all things:
  const SphereList& get_spheres() const{
    return spheres;
  }

 const LightSources& get_lights() const{
   return lightSource_list;
  }

 const ModelObjects& get_models() const{
   return modelObject_list;
 }
  
  // class instance variables:
 protected:
  
  Color background_color;
  
  // ambient illumination in the scene:
  string ambient_header;
  Color ambient_color;

  // Actual Lists:
  ModelObjects modelObject_list;
  LightSources lightSource_list;
  SphereList spheres;
  
};

#endif // SCENE_H_INCLUDE
