// PA4 Assignment
// Author: Tyrus Malmstrom
// Date  : 11/1/2016
// Class : CS410
// Email : tyrus.alexander.malmstrom@gmail.com
// 'Lights Up'

/*

                                       ~QUICK NOTE~

In building my raytracer for this assignment, I used a third party library called Eigen.

As taken from their official webpage, 
"Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms."

Furthermore, Here is a general overview of Eigen:

                                         ~Overview~
- Eigen is versatile.
  -It supports all matrix sizes, from small fixed-size matrices to arbitrarily large dense matrices, and even sparse matrices.
  -It supports all standard numeric types, including std::complex, integers, and is easily extensible to custom numeric types.
  -It supports various matrix decompositions and geometry features.
  -Its ecosystem of unsupported modules provides many specialized features such as non-linear optimization, matrix functions, a polynomial solver, FFT, and much more.

- Eigen is fast.
  -Expression templates allow to intelligently remove temporaries and enable lazy evaluation, when that is appropriate.
  -Explicit vectorization is performed for SSE 2/3/4, ARM NEON (32-bit and 64-bit), PowerPC AltiVec/VSX (32-bit and 64-bit) instruction sets, and now S390x SIMD (ZVector) with graceful fallback to non-vectorized code.
  -Fixed-size matrices are fully optimized: dynamic memory allocation is avoided, and the loops are unrolled when that makes sense.
  -For large matrices, special attention is paid to cache-friendliness.

- Eigen is reliable.
  -Algorithms are carefully selected for reliability. Reliability trade-offs are clearly documented and extremely safe decompositions are available.
  -Eigen is thoroughly tested through its own test suite (over 500 executables), the standard BLAS test suite, and parts of the LAPACK test suite.

- Eigen is elegant.
  -The API is extremely clean and expressive while feeling natural to C++ programmers, thanks to expression templates.
  -Implementing an algorithm on top of Eigen feels like just copying pseudocode.
  -Eigen has good compiler support as we run our test suite against many compilers to guarantee reliability and work around any compiler bugs. Eigen also is standard C++98 and maintains very reasonable compilation times.

                                         ~Requirements~

  Eigen doesn't have any dependencies other than the C++ standard library.

  We use the CMake build system, but only to build the documentation and unit-tests, and to automate installation. If you just want to use Eigen, you can use the header files right away. There is no binary library to link to,
  and no configured header file. Eigen is a pure template library defined in the headers.


~More information can be found here: http://eigen.tuxfamily.org/index.php?title=Main_Page

***************************************************************************************************

I also used another third party library called Tiny obj Loader.

As taken from their official webpage, 
"Tiny but powerful single file wavefront obj loader written in C++. No dependency except for C++ STL. It can parse 10M over polygons with moderate memory and time."

~More information may be found here: https://syoyo.github.io/tinyobjloader/


Best,
Tyrus Malmstrom.
  
*/

// include headers:
#include <iostream>
#include <string>
#include "Camera.h"

using namespace std;

int main(int argc, char *argv[]){

  // Instaniating class object(s):
  int numSpheres,numModels;
  Camera camera;
  camera.parseScene( argv[1] ); // this function does a whole lot, look at the source code for more details.  
  camera.buildRM();
  camera.calculateRays();

  numSpheres = camera.numberOfSpheres();
  // cout << "main numSpheres = " << numSpheres << endl;
  numModels  = camera.numberOfModels();
  // cout << "main numModels = " << numModels << endl;

  // 'special' case for scene4:
  if( numSpheres == 2 && numModels == 1){
    // cout << "Sphere and Model " << endl;
    camera.rayTriangleIntersection(); // get model stuff
    camera.writeSpheresAndModels( argv[2] );
  }
  else if( numSpheres > 0 ){
    // cout << "Spheres" << endl;
    camera.writeSpheres( argv[2] );
  }
  else if( numModels > 0 ){
    // cout << "Models" << endl;
    camera.rayTriangleIntersection();  
    camera.writeModels( argv[2] );
  }

  
  return 0;
}
