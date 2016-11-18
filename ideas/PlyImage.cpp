// Tyrus Malmstrom
// 11/1/2016
// PlyImage.cpp class for reading the .ply 'header'


// directives:
#include <iostream>
#include <string>
#include <vector>
#include <fstream> // Stream class to both read and write from/to files.
#include <sstream>
#include <algorithm> // std::remove
#include <iterator> // for iterator
#include <Eigen/Dense>
#include "PlyImage.h"
#include "ModelObject.h"
#include "Face.h"

// namespace
using namespace std;
using Eigen::MatrixXd;
// using Eigen::RowVectorXd;
using Eigen::Vector4d;

// functin declarations:
int remove_comments(const string &current_line);
int get_properties(const string &current_line);
void print_comments(const vector<string> &comments);
void print_properties(const vector<string> &properties);

// Macros:
#define DEBUG false
#define COMMENT_KEYWORD "comment"
#define PROPERTY_KEYWORD "property"
#define END_HEADER "end_header"

void PlyImage::readPlyFile(const string& fileName, ModelObject& obj, Face& faces){

  string line;
  ifstream infile(fileName);

  // checking to see if the file opened:
  if( !infile.is_open() ){
    // cout << "Sorry! Could open " << fileName << "!" << endl;
  }

  // grab the file header:
  getline(infile, file_header);
  if(DEBUG) cout << "The file header is: " << file_header <<  endl;

  // grab the file format :: ascii or binary
  getline(infile, file_format);
  if(DEBUG) cout << "The file format is: " << file_format << endl;


  // now dealing with comments, there may be one, or there may be more:
  stringstream comment_stream;
  string first_word;
  number_of_comments = 0;
  do{
    getline(infile, line);
    comment_stream << line;
    comment_stream >> first_word;
    // if(first_word  == COMMENT_KEYWORD) cout << "comment (" << number_of_comments << ") is: " << line << endl;
    if(first_word  == COMMENT_KEYWORD) number_of_comments++;
    if(first_word  == COMMENT_KEYWORD) comment_vector.push_back(line);
    comment_stream.str( string() );
    comment_stream.clear(); // clearing the stringstream
  }while( remove_comments(first_word) );

  // cout << "The number of comments in the .ply file are: " << number_of_comments << endl;
  // print_comments(comment_vector);

  stringstream element_line;
  string discard;
  // grab what type of element is in the ply file:
  element_line << line;
  first_element_def = element_line.str();
  // cout << "element_line.str() " << element_line.str() << endl;
  element_line >> discard; // this is just the word 'element'/
  element_line >> type_of_element;
  if(DEBUG) cout << "\"" << type_of_element << "\"" << " is the element defined in the .ply file" << endl;
  element_line >> number_of_verticies; // how many verticies
  if(DEBUG) cout << number_of_verticies << " verticies in the file." << endl;
  element_line.str( string() );
  element_line.clear(); // clearing the stringstream.

  obj.set_verticies(number_of_verticies);

  // cout << line << endl;


  // get the properties of the .ply file:
  stringstream property_stream;
  string firstt_word;
  number_of_properties = 0;
  do{
    getline(infile, line);
    property_stream << line;
    property_stream >> firstt_word;
    // cout << "firstt_word = " << firstt_word << endl;
    // if(firstt_word == PROPERTY_KEYWORD) cout << "property (" << number_of_properties << ") is: " << line << endl;
    if(firstt_word == PROPERTY_KEYWORD) number_of_properties++;
    if(firstt_word == PROPERTY_KEYWORD) properties_vector.push_back(line);
    property_stream.str( string() );
    property_stream.clear();
  }while(get_properties(firstt_word) );

  // cout << line << endl;

  string discardd;
  // grab what type of element it is:
  element_line << line;
  second_element_def = element_line.str();
  element_line >> discardd >> type_of_element2; // discarding the 'element' word and storing the element itself
  if(DEBUG) cout << "\"" << type_of_element2 << "\"" << " is another element defined in the .ply file" << endl;
  element_line >> number_of_faces; // how many faces (polygons) in the .ply file.
  if(DEBUG) cout << number_of_faces << " faces are in the file." << endl;

  // setting the ModelObject properties:
  obj.set_faces(number_of_faces);

  // cout << line << endl;

  do{
    getline(infile, line);
    property_stream << line;
    property_stream >> firstt_word;
    // cout << "firstt_word = " << firstt_word << endl;
    // if(firstt_word == PROPERTY_KEYWORD) cout << "property (" << number_of_properties << ") is: " << line << endl;
    if(firstt_word == PROPERTY_KEYWORD) number_of_properties++;
    if(firstt_word == PROPERTY_KEYWORD) properties_vector.push_back(line);
    property_stream.str( string() );
    property_stream.clear();
  }while(get_properties(firstt_word) );

  // cout << "The number of properties in the .ply file are: " << number_of_properties << endl;
  // print_properties(properties_vector);

  stringstream end;
  end << line;
  end >> end_header;
  // line at this point should be at the end of the hearder:
  if( end_header == END_HEADER ){
    // cout << "At the end of the header of the .ply file!" << endl;
  }


  // cout << "the line = " << line << endl;

  // start to read the data into the ModelObject obj:
  readData(infile, obj, faces);

}

void PlyImage::readData(ifstream& istr, ModelObject& obj, Face& faces){

  // reading in the verticies from the file:
  int row_verts = obj.get_verticies();

  // Initializing Matrix size:
  vertices.resize(row_verts,3);

  // cout << "rows_verts = " << rows_verts << endl;
  // reading in the list of verticies:
  for(int r = 0; r < row_verts; r++){
    for(int c = 0; c < 3; c++){
      if( !(istr >> vertices(r,c)) ){
	// cout << "Couldn't read in all the vertices! " << endl;
      }
    }
  }

  // cout << vertices << "\n" << endl;


  // set ModelObject vertex list -->  VERY IMPORTANT:
  obj.set_vertex_list( vertices );

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // copying what's left in the ifstream into a stringstreamm:
  /* OLD WAY
  stringstream sout;
  copy(istreambuf_iterator<char>(istr),
       istreambuf_iterator<char>(),
       ostreambuf_iterator<char>(sout));

  all_faces = sout.str();
  // cout << all_faces << endl;
  */
  
  int cols = obj.get_faces();
  list_faces.resize(4,cols);
  for(int i = 0; i < cols; i++){
    for(int c = 0; c < 4; c++){
      istr >> list_faces(c,i);
    }
  }

  // cout << list_faces << endl;
  obj.set_faces_list( list_faces );

  // allocate space for Face vector:
  F = vector< Face >( list_faces.cols() );
  for(int i = 0; i < list_faces.cols(); i++){
    Vector4d tmpface = list_faces.col(i);
    // cout << tmpface << endl;
    F[i] = Face( tmpface(1), tmpface(2), tmpface(3) );
    F[i].map(vertices);
    faces.addFace(F[i]);
    // F[i].pprint();
  }

  
  /*while( getline(istr,line) ){
    fn << line;
    fn >> fnn;
    RowVectorXd current_row( fnn );
    for(int i = 0; i < fnn; i++){
      fn >> index;
      current_row(i) = index;
    } // end for.
    cout << current_row << endl;
    // clearing stringstream:
    fn.str( string() );
    fn.clear();
  } // end while.
  */

  // another way I wanted to read in the triangle faces:
  /*int face_cols = 4;
  int face_rows = obj.get_faces();
  list_faces.resize(face_rows,face_cols);

  for(int i = 0; i < list_faces.rows(); i++){
    for(int j = 0; j < list_faces.cols(); j++){
      if( !( sout >> list_faces(i,j) ) ){
      	// cout << "Failed read when reading in the list of verticies" << endl;
      }
    }
  }

  obj.set_faces_list( list_faces );
  */

}

// writing out to the disk with the file name:
void PlyImage:: writePlyFile(ModelObject& obj, const string& fileName, const string& fileSuffix){

  // This just grabs the 'name' of the file before the .ply extension.
  string name_bf_file_exstension = fileName.substr( 0, (fileName.length() - 4) );

  ofstream outfile( name_bf_file_exstension + fileSuffix );

  if( !outfile.is_open() ){
    // cout << "Sorry! Couldn't write out to the file." << endl;
  }

  // start writing out to the file:
  outfile << file_header << endl << file_format << endl;
  // comments:
  for(int i = 0; i < static_cast<int>(comment_vector.size()); i++){
    outfile << comment_vector[i] << endl;
  }
  outfile << first_element_def << endl;
  // properties
  for(int i = 0; i < static_cast<int>(properties_vector.size()-1 ); i++){
    outfile << properties_vector[i] << endl;
  }
  outfile << second_element_def << endl << properties_vector[properties_vector.size()-1] << endl << end_header << endl;
  // write the verticies:
  MatrixXd write_out = obj.get_main_vertex_list();
  for(int i = 0; i < write_out.rows(); i++){
    for(int j = 0; j < write_out.cols(); j++){ // should be three columns.
      if( j == 2 ){
	outfile << write_out(i,j);
      }else{
	outfile << write_out(i,j) << " ";
      }
    }
    if( i != (write_out.cols()-1) ){
      outfile << endl;
    }
  }

  outfile << all_faces;
  outfile.close();

}

int remove_comments(const string &current_line){
  // cout << "In the remove_comments function and the current line is: " << current_line << endl;
  if( current_line == "comment"){ // comment is the first word in the string.
    return 1;
  }
  return 0;
}

int get_properties(const string &current_line){
  // cout << "In the get properties function and the current line is: " + current_line << endl;
  if( current_line == "property"){
    return 1;
  }
  return 0;
}


void print_comments(const vector<string> &container_of_comments){
  cout << "==============================" << endl;
  cout << "Here are the comments:" << endl;
  for(int i = 0; i < static_cast<int>( container_of_comments.size() ); i++){
    cout << container_of_comments[i] << endl;
  }
  cout << "==============================" << endl;
}

void print_properties(const vector<string> &container_of_properties){
  cout << "==============================" << endl;
  cout << "Here are the properties:" << endl;
  for(int i = 0; i < static_cast<int>( container_of_properties.size() ); i++){
    cout << container_of_properties[i] << endl;
  }
  cout << "==============================" << endl;
}
