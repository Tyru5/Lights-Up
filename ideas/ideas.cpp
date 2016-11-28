
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


// void Camera::getModelFaces() {

//   /*
//     This member function gets the point of intersection for each face before getting the color for each face.
//    */

//   for( int i = 0; i < width; i++){
//     for( int c = 0; c < height; c++){ // for each ray

//       for( int m = 0; m < static_cast<int>(modelObject_list.size()); m++){ // for all sphere in scene
// 	modelObject_list[m].rayTriangleIntersection( Rays[i][height - c -1] );
//       }

//     } 
//   } // end of rays.

//   cout << "~~~~Printing out the faces from the Camera class~~~~" << endl;
//   for( int m = 0; m < static_cast<int>(modelObject_list.size()); m++){ // for all sphere in scene
//     modelObject_list[m].printFaces();
//   }
  
// }

/*
void Camera::getModelFacesRGB() {

  for( int i = 0; i < width; i++){
    for( int c = 0; c < height; c++){ // for each ray
      
      for( int m = 0; m < static_cast<int>(modelObject_list.size()); m++){ // for all sphere in scene
	int number_of_faces = modelObject_list[m].numberOfFaces();
	// cout << "number of faces in loop= " << number_of_faces << endl;
	for(int f = 0; f < number_of_faces; f++){

	  tuple<bool, Color> bool_color = modelObject_list[m].getRayModelRGB( Rays[i][height - c -1], modelObject_list[m].getFace(f), ambient_color, lightSource_list );
	  // cout << "Color from getModelFacesRGB = " << color;
	  model_colors.push_back( bool_color );
	} // end of faces
	
      } // end of models
      
    } 
  } // end of rays.

  
  cout << "Printing the colors in the camera class" << endl;
  for ( const auto& i : model_colors ) {
    cout << get<0>(i) << get<1>(i);
  }

}
*/
