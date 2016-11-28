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

// Algorithm for Ray Triangle Intersection:
void Camera::computeDist( const Face& current_face ){

  /*For each pixel, throw ray out of focal point
    and calculate colour along ray;
    Fill in pixel value;
  */

  // int num_faces = obj.get_faces();
  // cout << "Polygon count: " << num_faces << endl;
  
  double beta;
  double gamma;
  double t;
  
  Vector3d origin(0,0,0);
  Vector3d direction(0,0,0); // origin, direction
  
  Vector3d A(0,0,0);
  Vector3d B(0,0,0);
  Vector3d C(0,0,0); // face vertices

  Matrix3d mtm(3,3);
  Matrix3d Mx1,Mx2,Mx3;
  double detMTM, detMTM1, detMTM2, detMTM3;
  
 

  // cout << faces << endl;
  
  for(int i = 0; i < width; i++){ // for each pixel on the image plane...
    for(int c = 0; c < height; c++){
      
      origin = Rays[i][c].origin;
      // cout << "O = \n" << O << endl;
      direction = Rays[i][c].direction;
      // cout << "D = \n" << D << endl;
      
     
      A = current_face.getA();
      // cout << "A = \n" << A << endl;
      B = current_face.getB();
      // cout << "B = \n" << B << endl;
      C = current_face.getC();
      // cout << "C = \n" << C << endl;
      
      // Find vectors for two edges sharing V1 (which is A in my case):
      Vector3d AB = A-B;
      Vector3d AC = A-C;
      Vector3d al = A-origin;
      
      mtm.col(0) = AB;
      mtm.col(1) = AC;
      mtm.col(2) = direction;
      
      // cout << mtm << endl;
      
      detMTM = mtm.determinant();
      
      Mx1 = mtm;
      Mx2 = mtm;
      Mx3 = mtm;
      
      Mx1.col(0) = al;  
      detMTM1 = Mx1.determinant();
      
      Mx2.col(1) = al;
      detMTM2 = Mx2.determinant();
      
      Mx3.col(2) = al;
      detMTM3 = Mx3.determinant();
      
      beta  = detMTM1/detMTM;
      // cout << "Beta: " << beta << endl;      
      gamma = detMTM2/detMTM;
      // cout << "Gamma: " << gamma << endl;
      t     = detMTM3/detMTM;
      // cout << " computed t: = " << t << endl;
      
      // Error Checking:
      if( beta >= 0.0 && gamma >= 0.0 && (beta+gamma <= 1.0) && t >= 0.0){ // ray intersect!
	// cout << "Ray intersected with face!" << endl;
	// cout << " computed t intersected: = " << t << endl;
	// cout << "Beta: " << beta << endl;
	// cout << "Gamma: " << gamma << endl;
	
	// checking t val:
	if( t <= ts[i][c] || ts[i][c] == -1.0 ){
	  ts[i][c] = t; // smallest tval
	  ptof[i][c] = origin + t* direction; // smallest point of intersection on face.
	}
	
      }
      
    }// end inner for loop.
  }// end outer for loop.

}


void Camera::rayTriangleIntersection(){

  int number_of_faces = modelObject_list[0].numberOfFaces();
  // allocate space for ts:
  ts = vector< vector< double > >(width, vector<double>( height, -1.0 )  );
  ptof= vector< vector< Vector3d > >(width, vector<Vector3d>( height, Vector3d(-1,-1,-1) )  );

  // print_ts(ts);
  // cout << "before" << endl;
  // print_ptof();
  
  for(int i = 0; i < number_of_faces; i++){
    // cout << face.getFace(i) << endl;
    computeDist( modelObject_list[0].getFace(i) ); // for now, only one model so 0th model in the vector of models
  }

  // print_ts(ts);
  // cout << "after " << endl;
  // print_ptof();
  // cout << "Polygon count (faces) : " << number_of_faces << endl;
  
}


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


void Camera::writeSpheresAndModels( const string& out_file ){

  ofstream out( out_file );
  if( !out ) cout << "Sorry! Couldn't write out the file: " << out_file << endl;

  sphere_model_pixs = vector< vector<RowVector3i> >(width, vector<RowVector3i>(height, RowVector3i(0,0,0) ) ); // pretty awesome
  // printPixs();
  
  // start writing out to the file:
  out << "P3 " << endl;
  out << width << " " << height << " 255" << endl;

  // Map pixel, get only one *very important :: talked with jake*
  for( int i = 0; i < width; i++){
    for( int c = 0; c < height; c++){ // for each ray

      // Spheres first:
      for(int sp = 0; sp < static_cast<int>(spheres.size()); sp++){
	
	tuple<bool, Color> res = spheres[sp].getRaySphereRGB( Rays[i][height - c -1], ambient_color, lightSource_list );
	if( get<0>(res) ){
	  sphere_model_pixs[i][c] = mapColour( get<1>(res) );
	  // cout << "Sphere pix[i,c] = " << sphere_model_pixs[i][c] << endl;
	}
      }

      // Now models:      
      for( int m = 0; m < static_cast<int>(modelObject_list.size()); m++){ // for all sphere in scene
      	int number_of_faces = modelObject_list[m].numberOfFaces();
      	// cout << "number of faces in loop= " << number_of_faces << endl;
      	for(int f = 0; f < number_of_faces; f++){
	  
      	  tuple<bool, Color> res  = modelObject_list[m].getRayModelRGB( Rays[i][height - c -1], modelObject_list[m].getFace(f), ptof[i][height - c - 1], ambient_color, lightSource_list );

      	  if( get<0>(res) ){
      	    sphere_model_pixs[i][c] = mapColour( get<1>(res) );
      	    // cout << "modelpix[i,c] = "<< i << ", " << c << ", " << sphere_model_pixs[i][c] << endl;
      	  }
   
      	} // end of model faces
      } // end of models
      
    } 
  } // end of rays.

  // printPixs();
  
  // now writing out:
  for(int i = 0; i < width; i++, out << endl){
    for(int j = 0; j < height; j++){ 
      // cout << model_pixs[j][i] << " "; // reversed to print out correctly
      out << sphere_model_pixs[j][i] << " "; // reversed to print out correctly
    }
  }
  
  out.close();

}


void Camera::writeSpheres( const string& out_file ){

  ofstream out( out_file );
  if( !out ) cout << "Sorry! Couldn't write out the file: " << out_file << endl;

  sphere_pixs = vector< vector<RowVector3i> >(width, vector<RowVector3i>(height, RowVector3i(0,0,0) ) ); // pretty awesome
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
	  sphere_pixs[i][c] = mapColour( get<1>(res) );
	  // cout << "pix[i,c] = " << pixs[i][c] << endl;
	  
	}
	
      }
    }
  }

  // now writing out:
  for(int i = 0; i < width; i++, out << endl){
    for(int j = 0; j < height; j++){
      out << sphere_pixs[j][i] << " "; // reversed to print out correctly
    }
  }
  
  out.close();

}


void Camera::writeModels( const string& out_file ){

  ofstream out( out_file );
  if( !out ) cout << "Sorry! Couldn't write out the file: " << out_file << endl;

  model_pixs = vector< vector<RowVector3i> >(width, vector<RowVector3i>(height, RowVector3i(0,0,0) ) ); // pretty awesome
  // printPixs();
  
  // start writing out to the file:
  out << "P3 " << endl;
  out << width << " " << height << " 255" << endl;

  // Map pixel, get only one *very important :: talked with jake*
  for( int i = 0; i < width; i++){
    for( int c = 0; c < height; c++){ // for each ray
      
      for( int m = 0; m < static_cast<int>(modelObject_list.size()); m++){ // for all sphere in scene
	int number_of_faces = modelObject_list[m].numberOfFaces();
	// cout << "number of faces in loop= " << number_of_faces << endl;
	for(int f = 0; f < number_of_faces; f++){
	  
	  tuple<bool, Color> res  = modelObject_list[m].getRayModelRGB( Rays[i][height - c -1], modelObject_list[m].getFace(f), ptof[i][height - c - 1], ambient_color, lightSource_list );

	  if( get<0>(res) ){
	    model_pixs[i][c] = mapColour( get<1>(res) );
	    // cout << "pix[i,c] = "<< i << ", " << c << ", " << model_pixs[i][c] << endl;
	  }
   
	} // end of faces
      } // end of models
      
    } 
  } // end of rays.

  // printPixs();
  
  // now writing out:
  for(int i = 0; i < width; i++, out << endl){
    for(int j = 0; j < height; j++){ 
      // cout << model_pixs[j][i] << " "; // reversed to print out correctly
      out << model_pixs[j][i] << " "; // reversed to print out correctly
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
      cout << model_pixs[i][c] << " ";
    }
    cout << endl;
  }

}

void Camera::print_ptof(){
  for(int i =  0; i < width; i++){
    for(int c = 0; c < height; c++){
      cout << ptof[i][c].transpose() << " ";
    }
    cout << endl;
  }
}
