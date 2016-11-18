// Tyrus Malmstrom
// Header file for the PlyImage.cpp

#ifndef PLYIMAGE_H_INCLUDE
#define PLYIMAGE_H_INCLUDE

// directives:
#include <iostream>
#include <string>
#include <vector>
#include "ModelObject.h"
#include <Eigen/Dense>
#include "Face.h"

// namespace:
using Eigen::MatrixXd;

class PlyImage{

 public:
  // constructor:
  PlyImage(){};
  // member functions:
  void readPlyFile(const std::string& fileName, ModelObject &obj, Face& faces);
  void readData(std::ifstream &ifstr, ModelObject &obj, Face& faces);
  void writePlyFile(ModelObject& obj, const std::string& fileName, const std::string& fileSuffix);



  // class instance variables:
 private:
  std::string file_header;
  std::string file_format;
  std::string first_element_def;
  std::string second_element_def;
  std::string type_of_element;
  std::string type_of_element2;
  std::string end_header;

  int number_of_comments;
  std::vector<std::string> comment_vector;

  int number_of_properties;
  std::vector<std::string> properties_vector;

  int number_of_verticies;
  MatrixXd vertices;

  int number_of_faces;
  MatrixXd list_faces;

  std::string all_faces;

  // to hold all the Face objects:
  std::vector< Face > F;


};

#endif // PLYIMAGE_H_INCLUDE
