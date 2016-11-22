// PA4 Assignment
// Author: Tyrus Malmstrom
// Date  : 11/1/2016
// Class : CS410
// Email : tyrus.alexander.malmstrom@gmail.com

/*

                                       ~QUICK NOTE~

In building my raytracer for this assignment, I used a third party library called Tiny obj Loader.

As taken from their official webpage, 
"Tiny but powerful single file wavefront obj loader written in C++. No dependency except for C++ STL. It can parse 10M over polygons with moderate memory and time."

More information may be found here: https://syoyo.github.io/tinyobjloader/

*/


// including directives:
#include <math.h> // for sqrt function
#include <algorithm> // replace

// Use this in *one* .cc
#define TINYOBJLOADER_IMPLEMENTATION // < -- *NEED to have the directives / includes like this*
#include "ModelObject.h"

// namespace
using namespace std;


// Macros:
#define DEBUG false


void ModelObject::pprint(ostream& out) const{
  out << "Model: " << obj_file << endl;
}

ostream& operator<< (ostream& out, const ModelObject& m){
  m.pprint( out );
  return out;
}


// Member function to parse the .obj model file:
void ModelObject::parseObj(){

  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, obj_file.c_str() );

  if (!err.empty()) { // `err` may contain warning message.
    cerr << err << endl;
  }
  
  if (!ret) {
    exit(1);
  }

 
}

/* ~Important note~
   I used this function and modified it for my own needs from syoyo @ https://github.com/syoyo/tinyobjloader
   - syoyo matains the project: Tiny Obj Loader, Tiny but powerful singile file wavefront obj loader program. 

*/
void ModelObject::PrintInfo() const{
  cout << "# of vertices  : " << (attrib.vertices.size() / 3) << endl;
  cout << "# of normals   : " << (attrib.normals.size() / 3)  << endl;
  cout << "# of texcoords : " << (attrib.texcoords.size() / 2) <<endl;

  cout << "# of shapes    : " << shapes.size() << endl;
  cout << "# of materials : " << materials.size() << endl;

  for (size_t v = 0; v < attrib.vertices.size() / 3; v++) {
    printf("  v[%ld] = (%f, %f, %f)\n", static_cast<long>(v),
           static_cast<const double>(attrib.vertices[3 * v + 0]),
           static_cast<const double>(attrib.vertices[3 * v + 1]),
           static_cast<const double>(attrib.vertices[3 * v + 2]));
  }

  for (size_t v = 0; v < attrib.normals.size() / 3; v++) {
    printf("  n[%ld] = (%f, %f, %f)\n", static_cast<long>(v),
           static_cast<const double>(attrib.normals[3 * v + 0]),
           static_cast<const double>(attrib.normals[3 * v + 1]),
           static_cast<const double>(attrib.normals[3 * v + 2]));
  }

  for (size_t v = 0; v < attrib.texcoords.size() / 2; v++) {
    printf("  uv[%ld] = (%f, %f)\n", static_cast<long>(v),
           static_cast<const double>(attrib.texcoords[2 * v + 0]),
           static_cast<const double>(attrib.texcoords[2 * v + 1]));
  }

  // For each shape
  for (size_t i = 0; i < shapes.size(); i++) {
    printf("shape[%ld].name = %s\n", static_cast<long>(i),
           shapes[i].name.c_str());
    printf("Size of shape[%ld].indices: %lu\n", static_cast<long>(i),
           static_cast<unsigned long>(shapes[i].mesh.indices.size()));

    size_t index_offset = 0;

    assert(shapes[i].mesh.num_face_vertices.size() == shapes[i].mesh.material_ids.size());

    printf("shape[%ld].num_faces: %lu\n", static_cast<long>(i),
           static_cast<unsigned long>(shapes[i].mesh.num_face_vertices.size()));

    // For each face
    for (size_t f = 0; f < shapes[i].mesh.num_face_vertices.size(); f++) {
      size_t fnum = shapes[i].mesh.num_face_vertices[f];

      printf("  face[%ld].fnum = %ld\n", static_cast<long>(f),
             static_cast<unsigned long>(fnum));

      // For each vertex in the face
      for (size_t v = 0; v < fnum; v++) {
        tinyobj::index_t idx = shapes[i].mesh.indices[index_offset + v];
        printf("    face[%ld].v[%ld].idx = %d/%d/%d\n", static_cast<long>(f),
               static_cast<long>(v), idx.vertex_index, idx.normal_index,
               idx.texcoord_index);
      }

      printf("  face[%ld].material_id = %d\n", static_cast<long>(f),
             shapes[i].mesh.material_ids[f]);

      index_offset += fnum;
    }

    printf("shape[%ld].num_tags: %lu\n", static_cast<long>(i),
           static_cast<unsigned long>(shapes[i].mesh.tags.size()));
    for (size_t t = 0; t < shapes[i].mesh.tags.size(); t++) {
      printf("  tag[%ld] = %s ", static_cast<long>(t),
             shapes[i].mesh.tags[t].name.c_str());
      printf(" ints: [");
      for (size_t j = 0; j < shapes[i].mesh.tags[t].intValues.size(); ++j) {
        printf("%ld", static_cast<long>(shapes[i].mesh.tags[t].intValues[j]));
        if (j < (shapes[i].mesh.tags[t].intValues.size() - 1)) {
          printf(", ");
        }
      }
      printf("]");

      printf(" floats: [");
      for (size_t j = 0; j < shapes[i].mesh.tags[t].floatValues.size(); ++j) {
        printf("%f", static_cast<const double>(
                         shapes[i].mesh.tags[t].floatValues[j]));
        if (j < (shapes[i].mesh.tags[t].floatValues.size() - 1)) {
          printf(", ");
        }
      }
      printf("]");

      printf(" strings: [");
      for (size_t j = 0; j < shapes[i].mesh.tags[t].stringValues.size(); ++j) {
        printf("%s", shapes[i].mesh.tags[t].stringValues[j].c_str());
        if (j < (shapes[i].mesh.tags[t].stringValues.size() - 1)) {
          printf(", ");
        }
      }
      printf("]");
      printf("\n");
    }
  }

  for (size_t i = 0; i < materials.size(); i++) {
    printf("material[%ld].name = %s\n", static_cast<long>(i),
           materials[i].name.c_str());
    printf("  material.Ka = (%f, %f ,%f)\n",
           static_cast<const double>(materials[i].ambient[0]),
           static_cast<const double>(materials[i].ambient[1]),
           static_cast<const double>(materials[i].ambient[2]));
    printf("  material.Kd = (%f, %f ,%f)\n",
           static_cast<const double>(materials[i].diffuse[0]),
           static_cast<const double>(materials[i].diffuse[1]),
           static_cast<const double>(materials[i].diffuse[2]));
    printf("  material.Ks = (%f, %f ,%f)\n",
           static_cast<const double>(materials[i].specular[0]),
           static_cast<const double>(materials[i].specular[1]),
           static_cast<const double>(materials[i].specular[2]));
    printf("  material.Tr = (%f, %f ,%f)\n",
           static_cast<const double>(materials[i].transmittance[0]),
           static_cast<const double>(materials[i].transmittance[1]),
           static_cast<const double>(materials[i].transmittance[2]));
    printf("  material.Ke = (%f, %f ,%f)\n",
           static_cast<const double>(materials[i].emission[0]),
           static_cast<const double>(materials[i].emission[1]),
           static_cast<const double>(materials[i].emission[2]));
    printf("  material.Ns = %f\n",
           static_cast<const double>(materials[i].shininess));
    printf("  material.Ni = %f\n", static_cast<const double>(materials[i].ior));
    printf("  material.dissolve = %f\n",
           static_cast<const double>(materials[i].dissolve));
    printf("  material.illum = %d\n", materials[i].illum);
    printf("  material.map_Ka = %s\n", materials[i].ambient_texname.c_str());
    printf("  material.map_Kd = %s\n", materials[i].diffuse_texname.c_str());
    printf("  material.map_Ks = %s\n", materials[i].specular_texname.c_str());
    printf("  material.map_Ns = %s\n",
           materials[i].specular_highlight_texname.c_str());
    printf("  material.map_bump = %s\n", materials[i].bump_texname.c_str());
    printf("    bump_multiplier = %f\n", static_cast<const double>(materials[i].bump_texopt.bump_multiplier));
    printf("  material.map_d = %s\n", materials[i].alpha_texname.c_str());
    printf("  material.disp = %s\n", materials[i].displacement_texname.c_str());
    printf("  <<PBR>>\n");
    printf("  material.Pr     = %f\n", static_cast<const double>(materials[i].roughness));
    printf("  material.Pm     = %f\n", static_cast<const double>(materials[i].metallic));
    printf("  material.Ps     = %f\n", static_cast<const double>(materials[i].sheen));
    printf("  material.Pc     = %f\n", static_cast<const double>(materials[i].clearcoat_thickness));
    printf("  material.Pcr    = %f\n", static_cast<const double>(materials[i].clearcoat_thickness));
    printf("  material.aniso  = %f\n", static_cast<const double>(materials[i].anisotropy));
    printf("  material.anisor = %f\n", static_cast<const double>(materials[i].anisotropy_rotation));
    printf("  material.map_Ke = %s\n", materials[i].emissive_texname.c_str());
    printf("  material.map_Pr = %s\n", materials[i].roughness_texname.c_str());
    printf("  material.map_Pm = %s\n", materials[i].metallic_texname.c_str());
    printf("  material.map_Ps = %s\n", materials[i].sheen_texname.c_str());
    printf("  material.norm   = %s\n", materials[i].normal_texname.c_str());
    std::map<std::string, std::string>::const_iterator it(
        materials[i].unknown_parameter.begin());
    std::map<std::string, std::string>::const_iterator itEnd(
        materials[i].unknown_parameter.end());

    for (; it != itEnd; it++) {
      printf("  material.%s = %s\n", it->first.c_str(), it->second.c_str());
    }
    printf("\n");
  }
}


void ModelObject::getVertices(){

  // Allocate right amount of space for verts,verts-norm,faces:
  // r x c
  vertices.resize( static_cast<int>( (attrib.vertices.size() / 3) ), 3);
  int row_count = 0;
  for(int v = 0; v < static_cast<int>( (attrib.vertices.size() / 3) ); v++){
    vertices(row_count,0) = attrib.vertices[3 * v + 0];
    vertices(row_count,1) = attrib.vertices[3 * v + 1];
    vertices(row_count,2) = attrib.vertices[3 * v + 2];      
    row_count++;
  }

  if(DEBUG) cout << "vertices = \n" << vertices<< endl;
  
}

void ModelObject::getFaces(){

  /*1st, get vertices*/
  getVertices();

  
  // Allocate right amount of space for verts,verts-norm,faces:
  F = vector< Face >( shapes[0].mesh.num_face_vertices.size() ); // zero b/c only one shape will ever be in it.
  vector<double> index_holder(3);

  double index;
  face_material.resize(3,3);
  size_t index_offset = 0;

  // For each face
  for (size_t f = 0; f < shapes[0].mesh.num_face_vertices.size(); f++) { // could be zero because only ever 1 shape
    int material_row_counter = 0;
    size_t fnum = shapes[0].mesh.num_face_vertices[f]; // could be zero because only ever 1 shape

    // For each vertex in the face
    for (size_t v = 0; v < fnum; v++) {
      tinyobj::index_t idx = shapes[0].mesh.indices[index_offset + v];
      index = idx.vertex_index;
      index_holder[v] = index;
      /*
      printf("    face[%ld].v[%ld].idx = %d/%d/%d\n", static_cast<long>(f),
  	     static_cast<long>(v), idx.vertex_index, idx.normal_index,
  	     idx.texcoord_index);
      */
    }
    
    if(DEBUG) cout << "material for face " << f << " = " << materials[ shapes[0].mesh.material_ids[f] ].diffuse[2] << endl;

    // cout << "reading ambient..." << endl;
    face_material(material_row_counter, 0) =  materials[ shapes[0].mesh.material_ids[f] ].ambient[0];
    face_material(material_row_counter, 1) =  materials[ shapes[0].mesh.material_ids[f] ].ambient[1];
    face_material(material_row_counter, 2) =  materials[ shapes[0].mesh.material_ids[f] ].ambient[2];
    material_row_counter++;    
    //cout << "reading diffuse..."<< endl;
    face_material(material_row_counter, 0) =  materials[ shapes[0].mesh.material_ids[f] ].diffuse[0];
    face_material(material_row_counter, 1) =  materials[ shapes[0].mesh.material_ids[f] ].diffuse[1];
    face_material(material_row_counter, 2) =  materials[ shapes[0].mesh.material_ids[f] ].diffuse[2];
    material_row_counter++;
    // cout << "reading specular..."<< endl;
    face_material(material_row_counter, 0) =  materials[ shapes[0].mesh.material_ids[f] ].specular[0];
    face_material(material_row_counter, 1) =  materials[ shapes[0].mesh.material_ids[f] ].specular[1];
    face_material(material_row_counter, 2) =  materials[ shapes[0].mesh.material_ids[f] ].specular[2];
    if(DEBUG) cout << "k matrix is = " << face_material << endl;

    F[f] = Face( index_holder[0], index_holder[1], index_holder[2], face_material );
    F[f].map( vertices, face_material );
    if(DEBUG) cout << F[f];
    
    
    index_offset += fnum;
  }


}

void ModelObject::rayTriangleIntersection( const int& width, const int& height ){ // Essentially res.

  /*
    2nd, find 't' (or how much to travel out the ray)
    For each face in the model, pass it and find 't'
  */
  for(int i = 0; i < shapes[0].mesh.num_face_vertices.size(); i++){ // 0 b/c only will ever be one shape
    computeDist( F[i] ); // pass each face from the model.
  }

}

// Algorithm for Ray Triangle Intersection:
void ModelObject::computeDist( const Ray& ray ){

  /*For each pixel, throw ray out of focal point
    and calculate colour along ray;
    Fill in pixel value;
  */

  double beta;
  double gamma;
  double t;
  
  Vector3d O(0,0,0);
  Vector3d D(0,0,0); // origin, direction
  
  Vector3d A(0,0,0);
  Vector3d B(0,0,0);
  Vector3d C(0,0,0); // face vertices

  Matrix3d mtm(3,3);
  Matrix3d Mx1,Mx2,Mx3;
  double detMTM, detMTM1, detMTM2, detMTM3;
  
 

  for(int i = 0; i < width; i++){ // for each pixel on the image plane...
    for(int c = 0; c < height; c++){
      
      O = Rays[i][c].origin;
      // cout << "O = \n" << O << endl;
      D = Rays[i][c].direction;
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
      Vector3d al = A-O;
      
      mtm.col(0) = AB;
      mtm.col(1) = AC;
      mtm.col(2) = D;
      
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

      // ADDED: Early break cases:
      if( beta > 0.0 || gamma < 0.0 ) break;
      if( beta+gamma > 1.0) break;
      if( t < 0.0 ) break;
      
      // Error Checking:
      if( beta >= 0.0 && gamma >= 0.0 && (beta+gamma <= 1.0) && t >= 0.0){ // ray intersect!
	// cout << "Ray intersected with face!" << endl;
	// cout << " computed t intersected: = " << t << endl;
	// cout << "Beta: " << beta << endl;
	// cout << "Gamma: " << gamma << endl;
	
	// checking t val:
	if( t <= ts[i][c] || ts[i][c] == -1.0){
	  ts[i][c] = t;
	}
	
      }
      
    }// end inner for loop.
  }// end outer for loop.

}


// void ModelObject::print_ts(const vector<vector<double>>& vect){
//   for(int i =  0; i < width; i++){
//     for(int c = 0; c < height; c++){
//       cout << vect[i][c] << " ";
//     }
//     cout << endl;
//   }
// }
