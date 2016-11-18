// Tyrus Malmstrom
// 11/9/2016
// Scene.cpp class for handling the scene


// directives:
#include <iostream>
#include <string>
#include <vector>
#include <fstream> // Stream class to both read and write from/to files.
#include <sstream>
#include <Eigen/Dense>
#include "Scene.h"

// namespace:
using namespace std;
using Eigen::Vector4d;

// Macros:
#define DEBUG true
#define AMBIENT_LIGHT  "ambient"
#define LIGHT_KEYWORD  "light"
#define SPHERE_KEYWORD "sphere"
#define MODEL_KEYWORD  "model"


void Scene::parseScene( const string& scene_file, Camera& camera ){

  ifstream current_scene( scene_file );
  if( !current_scene ){
    cout << "Sorry! Could open " << scene_file << "!" << endl;
  }

  // Parse and 'build' the Camera:
  camera.parseCameraSpecs( scene_file, current_scene );
  camera.buildRM();

  // Start of processing the scene itself:
  string line;
  double r,g,b;
  double x,y,z,w;


  stringstream iss;
  string identifier;
  
  /*For lightSources*/
  Color eX;

  /*For Sphere*/
  Color materialX;
  double radiusX;

  /*For models*/
  double omega, uv,vv,wv; // axis angle rotation
  string model_obj_file; // name of .obj file for model(s)

  
  while( getline( current_scene, line ) ){
    iss << line;
    cout << iss.str() << endl;
    iss >> identifier;
    if(DEBUG) cout << "identifier = " << identifier << endl;

    // Grabbing the ambient light (only one in our class) in the scene:
    if( identifier == AMBIENT_LIGHT ){
      iss >> r >> g >> b;
      ambient_color = Color(r,g,b);
      if(DEBUG) cout << ambient_color << "(ambient color)" << endl;
    }

    else if( identifier == LIGHT_KEYWORD){
      /*
	After specifying the amount of ambient light in the scene, zero or more light sources may be specified.
	The first four values given are the x, y, z and w coordinates of the light source in world coordinates.
	The last three values indicate the red, green and blue levels of the light source on a zero-one scale.
      */      
      iss >> x >> y >> z >> w >> r >> g >> b;
      Vector4d pX(x,y,z,w);
      eX = Color(r,g,b);
      LightSource lightX( pX, eX );
      addLightSources( lightX );
    }

    else if( identifier == SPHERE_KEYWORD){
      // Following the light sources come zero or more spheres:
      iss >> x >> y >> z >> radiusX >> r >> g >> b;
      Vector3d centerX( x,y,z );
      materialX = Color( r,g,b );
      Sphere sphereX( centerX, radiusX, materialX  );
      addSphere( sphereX );
    }
    
    else if( identifier == MODEL_KEYWORD){
      /* 
	 Finally, zero or more polygonal models may be specfied for inclusion in the scene.
	 Note the first seven values may indicate a Model to World transformation
	 first three are the x, y, z translation from model to world coordinates
	 next four specify an axis-angle rotation
	 The last argument is a string indicating the name of the file containing the 3D polygonal model in OBJ format.
      */
      iss >> x >> y >> z  >> omega >> uv >> vv >> wv >> model_obj_file;
      // maybe have a function here that would set up the model / parse the obj file;

    }

    // clearing the stringstream:
    iss.str( string() );
    iss.clear();
    
    
  } // end of while statement.

  if( !model_obj_file.empty() ) cout << "The model object file is = "  << model_obj_file << endl;
  
  
}
