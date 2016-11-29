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

Vector3d getSnrm(const Face& face, const Vector3d& toCamerVector);

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

void ModelObject::getVnertices(){

  // cout << "how many vns = " << (attrib.normals.size() / 3) << endl;
  vn = vector< Vector3d > (static_cast<int>( (attrib.normals.size() / 3) ) );
  Vector3d vnt(3);
  for (size_t v = 0; v < attrib.normals.size() / 3; v++) {
    /*
      printf("  n[%ld] = (%f, %f, %f)\n", static_cast<long>(v),
      static_cast<const double>(attrib.normals[3 * v + 0]),
      static_cast<const double>(attrib.normals[3 * v + 1]),
      static_cast<const double>(attrib.normals[3 * v + 2]));
    */
    vnt(0) = attrib.normals[3 * v + 0];
    vnt(1) = attrib.normals[3 * v + 1];
    vnt(2) = attrib.normals[3 * v + 2];
    // cout << "vn " << v << " =\n" << vnt << endl;
    vn[v] = vnt;
  }


}

void ModelObject::getFaces(){

  /*1st, get vertices*/
  getVertices();
  getVnertices();
  
  // Allocate right amount of space for verts,verts-norm,faces:
  F = vector< Face >( shapes[0].mesh.num_face_vertices.size() ); // zero b/c only one shape will ever be in it.

  vector<double> index_holder(3);
  vector<double> vn_index_holder(3);

  face_material.resize(3,3);
  size_t index_offset = 0;

  // For each face
  for (size_t f = 0; f < shapes[0].mesh.num_face_vertices.size(); f++) { // could be zero because only ever 1 shape
    int material_row_counter = 0;
    size_t fnum = shapes[0].mesh.num_face_vertices[f]; // could be zero because only ever 1 shape
    // cout << "fnum = " << fnum << endl;

    // For each vertex in the face
    for (size_t v = 0; v < fnum; v++) {
      tinyobj::index_t idx = shapes[0].mesh.indices[index_offset + v];
      // index = idx.vertex_index;
      index_holder[v] = idx.vertex_index;
      vn_index_holder[v] = idx.normal_index;
      /*
	printf("    face[%ld].v[%ld].idx = %d/%d/%d\n", static_cast<long>(f),
	static_cast<long>(v), idx.vertex_index, idx.normal_index,
	idx.texcoord_index);
      */
    }
    
    if(DEBUG) cout << "material for face " << f << " = " << materials[ shapes[0].mesh.material_ids[f] ].name<<endl;

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
    // if(DEBUG) cout << "k matrix is = " << face_material << endl;
    
    F[f] = Face( index_holder[0], index_holder[1], index_holder[2], face_material, vn[ vn_index_holder[0] ] ); // should matter what normal index...
    F[f].map( vertices, face_material, vn[ vn_index_holder[0] ] ); // index for that face
    if(DEBUG) cout << F[f];
    
    
    index_offset += fnum;
  }

}

/* Found out that i don't need to do this! */
// tuple<bool, Color> ModelObject::getRayModelRGB( const Ray& ray, const Face& face, const Vector3d& ptof, const Color& ambl, const vector<LightSource>& lights ){

//   /*
//     Given a certain ray-triangle (face) intersection, compute the RGB off the surface:
//   */
  
//   // Now, for each face that each model has, get the RGB of each face
//   Vector3d zero1_vector(-1,-1,-1);
//   double alpha = 16.0;
//   Color color; // to start off, a blank color;

//   Vector3d toCameraVector  = ray.origin - ptof; toCameraVector = toCameraVector / toCameraVector.norm();  
  
//   // if( face.ptos != zero_vector ){ // if the ray did intersect with the triangle (face)
//   if( ptof != zero1_vector ){
    
//     Vector3d snrm = face.surface_normal; // getSnrm( face, toCameraVector ); // unit length <-- updated this!
//     // if(DEBUG) cout << "the snrm on sphere is = " << snrm.transpose() << " with ptos = " << ptos.transpose() << endl;
//     // Initial condition of the ambient lighting of the scene:
//     Vector3d fa = face.material.row(0); // zero row will be ambient
//     Color face_ambient = Color( fa(0), fa(1), fa(2) );
//     // cout << face_ambient << endl;
//     color = ambl * face_ambient;
//     // cout << color;
      
//     // cout << lights.size() << endl;
//     for( int z = 0; z < static_cast<int>( lights.size() ); z++){
    
//       Vector3d lp( lights[z].position(0), lights[z].position(1), lights[z].position(2) );
//       // if(DEBUG) cout << "light position = " << lp.transpose() << endl;
    
//       Vector3d toL = lp - ptof; toL = toL/toL.norm(); // unit length
//       // cout << "toL = " << toL.transpose() << " with associated ptos = " << face.ptos.transpose() << endl;
    
//       if( snrm.dot( toL ) > 0.0 ){ // checking now light behind the object

//   	Vector3d fd = face.material.row(1); // zero row will be ambient
//   	Color face_diffuse = Color( fd(0), fd(1), fd(2) );
// 	color += face_diffuse * lights[z].energy * snrm.dot( toL );
//   	// cout << "color2 = " << color;
//   	Vector3d toC  = ray.origin - ptof; toC = toC / toC.norm();
//   	// cout << "toC = " << toC.transpose() << " with associated ptos = " << ptos.transpose() << endl;
	
//   	Vector3d spR  = (2 * snrm.dot( toL ) * snrm) - toL;
//   	// cout << "spR = " << spR.transpose() << " with ptos of = " << ptos.transpose() << endl;

//   	// cout << toC.dot( spR ) << " ptos associated = " << ptos.transpose() << endl;; //<-- why not 16?

//   	Vector3d fs = face.material.row(2);
//   	Color face_specular = Color( fs(0), fs(1), fs(2) );
//   	color += face_specular * lights[z].energy *  pow( toC.dot( spR ), alpha );
//   	// cout << "color3 = " << color << "with ptos of = " << ptof.transpose() << endl;
	
//       }

//     }
    
//     // cout << "about to return the color." << endl;
//     return make_tuple(true, color);
    
//   }else{
    
//     return make_tuple(false, Color() );
      
//   }
    
// }


void ModelObject::printFaces() const{
  for(int i = 0; i < static_cast<int>( shapes[0].mesh.num_face_vertices.size() ); i++){
    cout << F[i] << endl;
  }
}

Vector3d getSnrm(const Face& face, const Vector3d& toCameraVector){

  Vector3d e1,e2,A,B,C,fnrm;
  A = face.getA();
  B = face.getB();
  C = face.getC();
  
  e1 = A-B;
  e2 = A-C;
  
  fnrm = e1.cross(e2);
  fnrm = fnrm / fnrm.norm();

  // condition check:
  if( fnrm.dot(toCameraVector) < 0.0 ) fnrm = -1 * fnrm; // flip the sign
  
  return fnrm;
}
