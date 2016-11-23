// Tyrus Malmstrom
// 11/1/2016
// Camera.cpp class for handling the camera specs.


// directives:
#include <fstream> // Stream class to both read and write from/to files.
#include <sstream>
#include <algorithm> // std::max
#include "Camera.h"

// namespace:
using namespace std;
using Eigen::RowVector3i;
using Eigen::Vector4d;
using Eigen::Vector3d;
using Eigen::Vector3i;


// Macros:
#define DEBUG false
#define AMBIENT_LIGHT  "ambient"
#define LIGHT_KEYWORD  "light"
#define SPHERE_KEYWORD "sphere"
#define MODEL_KEYWORD  "model"

void Camera::parseScene( const string& scene_file ){

  string line;
  double x,y,z,w;

  ifstream scene( scene_file );
  
  if( !scene ){
    cout << "Sorry! Could open " << scene_file << "!" << endl;
  }


  // Grab the first line:
  stringstream eye_stream;
  getline(scene, line);
  if(DEBUG) cout << "The first line in the camera file is: " << line << endl;

  eye_stream << line;
  eye_stream >> eye_header >> x >> y >> z;
  Vector3d a(x,y,z);
  EYE = a;
  if(DEBUG) cout << "The eye / focal point is at: \n" << EYE << endl;

  // grab lookap:
  stringstream lookap_stream;
  getline(scene, line);
  lookap_stream << line;
  lookap_stream >> look_header >> x >> y >> z;
  Vector3d b(x,y,z);
  LOOKAP = b;
  if(DEBUG) cout << "The look at point is at: \n" << LOOKAP << endl;


  // grab UPV:
  stringstream upv_stream;
  getline(scene, line);
  upv_stream << line;
  upv_stream >> supv >> x >> y >> z;
  Vector3d c(x,y,z);
  UPV = c;
  if(DEBUG) cout << "The up vector is at: \n" << UPV << endl;

  // grab distance away from image plane:
  stringstream dist_stream;
  getline(scene, line);
  dist_stream << line;
  dist_stream >> dist_header >> dist;
  // have to NEGATE the d, because we are looking DOWN the negative z axis:
  dist = -dist;
  if(DEBUG) cout << "The distance away from the image plane is: " << dist << "\n" <<  endl;

  // grab the bounds:
  stringstream bounds_stream;
  getline(scene, line);
  bounds_stream << line;
  bounds_stream >> bounds_header >> left >> bottom >> right >> top; // CHANGED THIS! Now it works! was initially reading in bounds wrong.
  if(DEBUG) cout << "Bottom, left, top, right is: " << left << " " <<  bottom << " " <<  right << " " <<  top << endl;

  // grab the resolution:
  stringstream res_stream;
  getline(scene, line);
  res_stream << line;
  res_stream >> res_header >> width >> height;
  if(DEBUG) cout << "width and height is: " << width << " " <<  " " << height << endl;

  cout << "Target resolution: " << width << " by " << height << endl;


  // Parsing the actual scene file:
  stringstream iss;
  string identifier;
  double red,green,blue;
  
  /*For lightSources*/
  Color eX;

  /*For Sphere*/
  Color materialX;
  double radiusX;

  /*For models*/
  double tx,ty,tz; // transformation
  double wx,wy,wz, theta; // axis angle rotation
  string model_obj_file; // name of .obj file for model(s)

  
  while( getline( scene, line ) ){
    iss << line;
    // cout << iss.str() << endl;
    iss >> identifier;
    if(DEBUG) cout << "identifier = " << identifier << endl;

    // Grabbing the ambient light (only one in our class) in the scene:
    if( identifier == AMBIENT_LIGHT ){
      iss >> red >> green >> blue;
      ambient_color = Color(red,green,blue);
      if(DEBUG) cout << ambient_color << "(ambient color)" << endl;
    }

    else if( identifier == LIGHT_KEYWORD){
      /*
	After specifying the amount of ambient light in the scene, zero or more light sources may be specified.
	The first four values given are the x, y, z and w coordinates of the light source in world coordinates.
	The last three values indicate the red, green and blue levels of the light source on a zero-one scale.
      */      
      iss >> x >> y >> z >> w >> red >> green >> blue;
      Vector4d pX(x,y,z,w);
      eX = Color(red,green,blue);
      LightSource lightX( pX, eX );
      addLightSources( lightX );
    }

    else if( identifier == SPHERE_KEYWORD){
      // Following the light sources come zero or more spheres:
      iss >> x >> y >> z >> radiusX >> red >> green >> blue;
      Vector3d centerX( x,y,z );
      materialX = Color( red,green,blue );
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
      iss >> tx >> ty >> tz >> wx >> wy >> wz >> theta >> model_obj_file;
      ModelObject modelX( model_obj_file, tx, ty, tz, wx, wy, wz, theta );
      addModels( modelX );
      // maybe have a function here that would set up the model / parse the obj file;

    }

    // clearing the stringstream:
    iss.str( string() );
    iss.clear();
    
    
  } // end of while statement.

  // cout << "How many \"models\" are in the file = " << modelObject_list.size() << endl;
  
  if( modelObject_list.size() != 0 ){

    // Now for each model, parse it and assign faces with cooresponding material props and vertex_norms;
    for( int i =0; i < static_cast<int>(modelObject_list.size()); i++){
      modelObject_list[i].parseObj();
      if(DEBUG) modelObject_list[i].PrintInfo();
      modelObject_list[i].getFaces();
      
    }
    
  } // end of if
  
}


void Camera::buildRM(){

  // Build Camera system origin and axes in world coordinates:
  /*
    Going to use the process described in Lecture Week 5:
    1) Point the z axis away --> Camera looks down the negative z axis:
    We have two points in 3D --> The eye and the look at point
    Gaze direction is L-E (however we are going to do E-L)
    So W axis of RM is going to be defined as: W = E-L/||E-L|| <-- make it unit length
  */

  WV = (EYE-LOOKAP);
  WV = WV/(WV.norm());
  if(DEBUG) cout << "W unit vector is: \n" << WV << endl;
  /* The U axis (horizontal axis) is perpendicular to a plane defined by UPV and W */
  UV = UPV.cross(WV);
  UV = UV/(UV.norm());
  if(DEBUG) cout << "U unit vector is: \n" << UV << endl;
  /*
    Given the first two axis, the third is:
    V = W X U
  */
  VV = WV.cross(UV);
  if(DEBUG) cout << "The V unit vector is: \n" << VV << endl;

}

// UPDATED:
void Camera::calculateRays(){

  Rays = vector< vector< Ray > >(width, vector<Ray>(height) );
  
  /*
    Code that creates a 3D point that represents a pixel on th image plane:
    As well as the rays. Get direction of each ray.
  */

  for(int i = 0; i < width; i++){
    for(int j = 0; j < height; j++){
      // cout << "i,j" << i << " " << j  << endl;
      
      double px = i/(width-1)  * (right-left) + left;
      double py = j/(height-1) * (top-bottom) + bottom;
      
      // Creating th pixel --> in world coordinates:
      // Awesome stuff man, vector + vector + vector + vector == point in the world.
      Vector3d pixelPoint = EYE + (dist * WV) + (px * UV) + (py * VV);
      // cout << "The pixel Point (3D point) in the world is: \n" << pixelPoint << endl;
      
      Vector3d rayd = pixelPoint - EYE;
      rayd = rayd/rayd.norm();
      Rays[i][j] =  Ray( pixelPoint, rayd );
      // Rays[i][j].pprint();
      
    }
  }

}

// void Camera::raySphereIntersection() {

//   for( int i = 0; i < width; i++){
//     for( int c = 0; c < height; c++){ // for each ray

//       for( int j = 0; j < static_cast<int>(spheres.size()); j++){ // for all sphere in scene

// 	tuple<bool, Color> bool_color = spheres[j].getRaySphereRGB( Rays[i][height - c -1], ambient_color, lightSource_list );
// 	// cout << "Color from raySphereIntersection = " << color;
// 	sphere_colors.push_back( bool_color );
	
//       }

//     } 
//   } // end of rays.


//   // cout << "Printing the colors in the camera class" << endl;
//   // for ( const auto& i : sphere_colors ) {
//   //   cout << get<0>(i) << get<1>(i);
//   // }
  
// }

RowVector3i Camera::mapColour( const Color &c ){

  // tuple(map(lambda(x) : ZZ(max(0,min(255,round(255.0 * x)))), res[1]));

  int red,green,blue;
  RowVector3i colorRGB(0,0,0);
  
  red   = max(0.0, min(255.0,round(255.0 * c.red )));
  green = max(0.0, min(255.0,round(255.0 * c.green)));
  blue  = max(0.0, min(255.0,round(255.0 * c.blue )));
  
  colorRGB(0) = red;
  colorRGB(1) = green;
  colorRGB(2) = blue;
  
  // cout << colorRGB << endl;
  return colorRGB;

}


void Camera::writeImage( const string& out_file ){

  ofstream out( out_file );
  if( !out ) cout << "Sorry! Couldn't write out the file: " << out_file << endl;

  pixs = vector< vector<RowVector3i> >(width, vector<RowVector3i>(height, RowVector3i(0,0,0) ) ); // pretty awesome
  // printPixs();
  
  // start writing out to the file:
  out << "P3 " << endl;
  out << width << " " << height << " 255" << endl;

  // Map pixel, get only one *very important :: talked with jake*
  Vector3i rgb(3);
  for(int i = 0; i < width; i++ ){
    for(int c = 0; c < height; c++ ){
      for(int sp = 0; sp < static_cast<int>(spheres.size()); sp++){

	tuple<bool, Color> res = spheres[sp].getRaySphereRGB( Rays[i][height - c -1], ambient_color, lightSource_list );
	if( get<0>(res) ){
	  // rgb = mapColour( sphere_colors[i * height + c] );
	  //rgb = mapColour( get<1>(res) );
	  // out << rgb(0) << " " << rgb(1) << " " << rgb(2) << " ";
	  pixs[i][c] = mapColour( get<1>(res) );
	  // cout << "pix[i,c] = " << pixs[i][c] << endl;
	  
	}
	
      }
    }
  }

  // now writing out:
  for(int i = 0; i < width; i++, out << endl){
    for(int j = 0; j < height; j++){
      out << pixs[j][i] << " "; // reversed to print out correctly
    }
  }
  
  out.close();

}


void Camera::writeImage2( const string& out_file ){

  ofstream out( out_file );
  if( !out ) cout << "Sorry! Couldn't write out the file: " << out_file << endl;

  pixs = vector< vector<RowVector3i> >(width, vector<RowVector3i>(height, RowVector3i(0,0,0) ) ); // pretty awesome
  // printPixs();
  
  // start writing out to the file:
  out << "P3 " << endl;
  out << width << " " << height << " 255" << endl;

  // Map pixel, get only one *very important :: talked with jake*
  Vector3i rgb(3);
  for(int i = 0; i < width; i++ ){
    for(int c = 0; c < height; c++ ){
      for(int m = 0; m < static_cast<int>(modelObject_list.size()); m++){

	tuple<bool, Color> res = modelObject_list[m].getRayModelRGB( width, height, Rays[i][height - c -1], ambient_color, lightSource_list );
	if( get<0>(res) ){
	  // rgb = mapColour( sphere_colors[i * height + c] );
	  //rgb = mapColour( get<1>(res) );
	  // out << rgb(0) << " " << rgb(1) << " " << rgb(2) << " ";
	  pixs[i][c] = mapColour( get<1>(res) );
	  // cout << "pix[i,c] = " << pixs[i][c] << endl;
	  
	}
	
      }
    }
  }

  // now writing out:
  for(int i = 0; i < width; i++, out << endl){
    for(int j = 0; j < height; j++){
      out << pixs[j][i] << " "; // reversed to print out correctly
    }
  }
  
  out.close();

}


// ==================HELPER FUNCTIONS=========================
void Camera::find_tmin_tmax(std::vector<std::vector<double>>& tvals){

  for(int i = 0; i < width; i++){
    for(int c = 0; c < height; c++){

      if(tvals[i][c] >= 0){

	if(tvals[i][c] < tmin) tmin = tvals[i][c];
	if(tvals[i][c] > tmax) tmax = tvals[i][c];
	
      }
      
    }
  }
  
}

// pa4
void Camera::printPixs() const{
  
  for(int i = 0; i < width; i++){
    for(int c = 0; c < height; c++){
      cout << pixs[i][c] << " ";
    }
    cout << endl;
  }

}
