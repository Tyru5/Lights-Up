// PA4 Assignment
// Author: Tyrus Malmstrom
// Date  : 11/1/2016
// Class : CS410
// Email : tyrus.alexander.malmstrom@gmail.com

// including directives:
#include <iostream>
#include <string>
#include <math.h> // for sqrt function
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm> // replace
#include <Eigen/Dense>
#include "ModelObject.h"

// namespace
using namespace std;
using Eigen::MatrixXd;

// Macros:
#define DEBUG true


void ModelObject::pprint(ostream& out) const{
  out << "Model: " << obj_file << endl;
}

ostream& operator<< (ostream& out, const ModelObject& m){
  m.pprint( out );
  return out;
}


// Member function to parse the .obj model file:
void ModelObject::parseObj(){

  string line;  
  ifstream obj( obj_file );
  if( !obj ) cout << "Sorry! Could open " << obj_file << "!" << endl;

  stringstream mdss;
  string identifier;
  // Parsing the whole obj file:
  while( getline( obj, line ) ){
    mdss << line;
    // cout << mdss.str() << endl;
    mdss >> identifier;
    // if(DEBUG) cout << "mdss identifier = " << identifier << endl;
    
    // Ignore comments, vt, and empty lines:
    if( identifier == "#" || identifier == "vt" || mdss.str().length() == 1 ){
      // clearing the stringstream:
      mdss.str( string() );
      mdss.clear();
      continue;
    }
    else if( identifier== "mtllib" ){
      /*This is the file that we have to load for the material properties of the face */
      cout << "found the mtllib file!" << endl;
      // have a function that does this:
    }
    else if( identifier == "v" ){
      number_of_vertices++;
    }
    else if( identifier == "vn" ){
      number_of_vertices_norm++;
    }
    else if( identifier == "f" ){
      number_of_faces++;
    }
    
    // clearing the stringstream:
    mdss.str( string() );
    mdss.clear();


  } // end of while loop
 

  // Allocate right amount of space for verts,verts-norm,faces:
  // r x c
  vertices.resize(number_of_vertices, 3);
  nertices.resize(number_of_vertices_norm,3);
  F= vector< Face >(number_of_faces);

  ifstream obj2( obj_file );

  /*For xyz vertices*/
  double x,y,z;
  double A,B,C;

  stringstream mdss2;
  string identifier2;
  
  int row_count = 0;
  int row_count2 = 0;
  int face_counter = 0;
  while( getline( obj2, line ) ){
    mdss2 << line;
    // cout << mdss2.str() << endl;
    mdss2 >> identifier2;
    // if(DEBUG) cout << "mdss identifier2 = " << identifier2 << endl;
    
    // Ignore comments, vt, and empty lines:
    if( identifier2 == "#" || identifier2 == "vt" || mdss2.str().length() == 1 ){
      // clearing the stringstream:
      mdss2.str( string() );
      mdss2.clear();
      continue;
    }

    else if( identifier2 == "v" ){
      mdss2 >> x >> y >> z;
      vertices(row_count,0) = x;
      vertices(row_count,1) = y;
      vertices(row_count,2) = z;      
      row_count++;
    }

    else if( identifier2 == "vn" ){
      mdss2 >> x >> y >> z;
      nertices(row_count2, 0) = x;
      nertices(row_count2, 1) = y;
      nertices(row_count2, 2) = z;
      row_count2++;
    }

    else if( identifier2 == "f" ){
      // cout << mdss2.str() << endl;      
      size_t offset = 0; // offset will be set to the length of characters of the "value" - 1.
      A = stod(&line[2], &offset); 
      B = stod(&line[offset + 6]); 
      C = stod(&line[offset + 12]);
      // cout << A << " " << B <<  " " << C << endl;
      F[face_counter] = Face( A-1, B-1, C-1 ); // had to subtract one for the offset
      face_counter++;
    }
  
    // clearing the stringstream:
    mdss2.str( string() );
    mdss2.clear();

    
  } // end of while loop

  if(DEBUG){
    cout << "vertices = \n" << vertices << endl;
    cout << "norms = \n" << nertices << endl;
  }
  
  // Map vertices to faces:
  for(int i = 0; i < static_cast<int>(F.size()); i++){
    F[i].map( vertices );
  }

  if(DEBUG){
    cout << "printing out the faces for the model:\n" << endl;
    for(int i = 0; i < static_cast<int>(F.size()); i++){
      cout << F[i] << endl;
    }
  }

 
}
