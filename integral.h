#ifndef INTEGRAL_H
#define INTEGRAL_H


//Macro d'acces rapide au voxel
#define integ(x,y,z) *(pin_tegral + (x) + (y)*int_egral_incs[1] + (z)*int_egral_incs[2] )
#define sourc(x,y,z) *(pso_urce   + (x) + (y)*sou_rce_incs[1]   + (z)*sou_rce_incs[2] )

#include <algorithm>  // req'd for std::min/max
#include <assert.h>
#include <vtkMultiThreader.h>

// undefine VS macros
#ifdef min
  #undef min
#endif

#ifdef max
  #undef max
#endif

#include <vtkImageData.h>



using namespace std;

typedef unsigned long long int ulli;

vtkImageData *ComputeIntegral(vtkImageData *img);
VTK_THREAD_RETURN_TYPE ThreadedIntegral (void *arg);






inline unsigned long long BoxIntegral(vtkImageData *img, int dim0, int dim1, int dim2, int size0, int size1, int size2)
{
  int* dims = img->GetDimensions();

  // The subtraction by one for row/col/lay is because row/col/lay is inclusive.
  int x1 = std::min(dim0, dims[0])  - 1;
  int y1 = std::min(dim1, dims[1])  - 1;
  int z1 = std::min(dim2, dims[2])  - 1;

  int x2 = std::min(dim0 + size0,   dims[0])  - 1;
  int y2 = std::min(dim1 + size1,   dims[1])  - 1;
  int z2 = std::min(dim2 + size2,   dims[2])  - 1;

  unsigned long long b222(0), b221(0), b212(0), b122(0), b112(0), b121(0), b211(0), b111(0);

  if (x2 >= 0 && y2 >= 0 && z2 >= 0)
		b222 = *static_cast<unsigned long long*>(img->GetScalarPointer(x2, y2, z2));
  if (x2 >= 0 && y2 >= 0 && z1 >= 0)
		b221 = *static_cast<unsigned long long*>(img->GetScalarPointer(x2, y2, z1));
  if (x2 >= 0 && y1 >= 0 && z2 >= 0)
		b212 = *static_cast<unsigned long long*>(img->GetScalarPointer(x2, y1, z2));
  if (x1 >= 0 && y2 >= 0 && z2 >= 0)
		b122 = *static_cast<unsigned long long*>(img->GetScalarPointer(x1, y2, z2));
  if (x1 >= 0 && y1 >= 0 && z2 >= 0)
		b112 = *static_cast<unsigned long long*>(img->GetScalarPointer(x1, y1, z2));
  if (x1 >= 0 && y2 >= 0 && z1 >= 0)
		b121 = *static_cast<unsigned long long*>(img->GetScalarPointer(x1, y2, z1));
  if (x2 >= 0 && y1 >= 0 && z1 >= 0)
		b211 = *static_cast<unsigned long long*>(img->GetScalarPointer(x2, y1, z1));
  if (x1 >= 0 && y1 >= 0 && z1 >= 0)
		b111 = *static_cast<unsigned long long*>(img->GetScalarPointer(x1, y1, z1));

  return std::max((unsigned long long)0, b222 - b221 - b212 - b122 + b112 + b121 + b211 - b111 );
}
//#define DEBUG


inline unsigned long long BoxIntegralOptim( const vtkImageData *img, const int& dim0, const int& dim1, const int& dim2, const int& size0, const int& size1, const int& size2, const long long int * int_egral_incs, int* dims,  const ulli* pin_tegral)
{

#ifdef DEBUG
	assert(dim0<dims[0]);
	assert(dim1<dims[1]);
	assert(dim2<dims[2]);

	assert(dim0 + size0<dims[0]);
	assert(dim1 + size1<dims[1]);
	assert(dim2 + size2<dims[2]);

  assert(dim0>0);
  assert(dim1>0);
  assert(dim2>0);

  assert(size0>0);
  assert(size1>0);
  assert(size2>0);
#endif

  // The subtraction by one for row/col/lay is because row/col/lay is inclusive.
  int x1 = dim0 - 1;
  int y1 = dim1 - 1;
  int z1 = dim2 - 1;

  int x2 = dim0 + size0 - 1;
  int y2 = dim1 + size1 - 1;
  int z2 = dim2 + size2 - 1;

#ifdef DEBUG
  assert(x2 >= 0 && y2 >= 0 && z2 >= 0);
  assert(x2 >= 0 && y2 >= 0 && z1 >= 0);
  assert(x2 >= 0 && y1 >= 0 && z2 >= 0);
  assert(x1 >= 0 && y2 >= 0 && z2 >= 0);
  assert(x1 >= 0 && y1 >= 0 && z2 >= 0);
  assert(x1 >= 0 && y2 >= 0 && z1 >= 0);
  assert(x2 >= 0 && y1 >= 0 && z1 >= 0);
  assert(x1 >= 0 && y1 >= 0 && z1 >= 0);
#endif

  return  integ(x2, y2, z2)
		- integ(x2, y2, z1)
		- integ(x2, y1, z2)
		- integ(x1, y2, z2)
		+ integ(x1, y1, z2)
		+ integ(x1, y2, z1)
		+ integ(x2, y1, z1)
		- integ(x1, y1, z1);
}



#endif
