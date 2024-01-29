/***********************************************************
*  --- OpenSURF ---                                       *
*  This library is distributed under the GNU GPL. Please   *
*  use the contact form at http://www.chrisevansdev.com    *
*  for more information.                                   *
*                                                          *
*  C. Evans, Research Into Robust Visual Features,         *
*  MSc University of Bristol, 2008.                        *
*                                                          *
************************************************************/

#include "surf.h"
#include <vtkMultiThreader.h>

#include <math.h>

using namespace std;

//-------------------------------------------------------
//! SURF priors (these need not be done at runtime)
const float pi = 3.14159f;

//-------------------------------------------------------

//! Constructor
Surf::Surf(vtkImageData *img, IpVec &ipts)
: ipts(ipts)
{
  this->img = img;

  dims = img->GetDimensions();
  incs = img->GetIncrements();
  pin_tegral = static_cast<ulli*>(img->GetScalarPointer());

}

//-------------------------------------------------------

//! Describe all features in the supplied vector
void Surf::getDescriptors( int radius, bool normalize )
{
  #pragma omp parallel for
  for (int i = 0; i < ipts.size(); ++i)
    this->getDescriptor(i, radius, normalize);

  return;
}


//! Describe all features in the supplied vector
void Surf::getRawDescriptors( int radius )
{
	#pragma omp parallel for
    for (int i = 0; i < ipts.size(); ++i)
      this->getRawDescriptor(i, radius);

	return;
}


//! Get the modified descriptor. See Agrawal ECCV 08
//! Modified descriptor contributed by Pablo Fernandez
void Surf::getDescriptor(int id, int radius, bool normalize)
{
  int y, x, z, sample_x, sample_y, sample_z, count=0;
  int i, j, k, xs, ys, zs;
  float ix, jx, kx;
  double scale, dx, dy, dz, mdx, mdy, mdz;
  double gauss_s1 = 0.f;
  double rx = 0.f, ry = 0.f, rz = 0.f, len = 0.f;

  Ipoint *ipt = &ipts[id];
  scale = ipt->scale;
  x = fRound(ipt->x);
  y = fRound(ipt->y);
  z = fRound(ipt->z);
  ipt->allocate(48);
  float *desc = ipt->descriptor;
  float halfRadius = (float ) ( radius - 1.0 ) / 2.0;

  i = -radius; //-10

  //Calculate descriptor for this key point
  while(i < radius) // 10
  {
    j = -radius; //-10

    while(j < radius) //10
    {
		k = -radius;

		while(k < radius)
		{

		  dx=dy=mdx=mdy=dz=mdz=0.f;

		  ix = ( float ) i + halfRadius; //Centre de chaque Bloc
		  jx = ( float ) j + halfRadius;
		  kx = ( float ) k + halfRadius;

		  xs = fRound( ipt->x + ( ix * scale ) );
		  ys = fRound( ipt->y + ( jx * scale ) );
		  zs = fRound( ipt->z + ( kx * scale ) );

		  for (int u = i; u < i + radius; ++u)
		  {
			for (int v = j; v < j + radius; ++v)
			{
				for (int w = k; w < k + radius; ++w)
				{
				  //Get coords of sample point on the rotated axis
				  sample_x = fRound(x + u*scale); // x + k
				  sample_y = fRound(y + v*scale); // y + l
				  sample_z = fRound(z + w*scale);

				  //Get the gaussian weighted x and y responses
				  gauss_s1 = gaussian( ( float ) xs-sample_x,
					 (float ) ys-sample_y,  (float ) zs-sample_z, 2.5f*scale);
          			//gauss_s1 = 1;
				  rx = gauss_s1*haarXOptim(sample_x, sample_y, sample_z, 2*fRound(scale));
				  ry = gauss_s1*haarYOptim(sample_x, sample_y, sample_z, 2*fRound(scale));
				  rz = gauss_s1*haarZOptim(sample_x, sample_y, sample_z, 2*fRound(scale));

				  dx += rx;
				  dy += ry;
				  dz += rz;
				  mdx += fabs(rx);
				  mdy += fabs(ry);
				  mdz += fabs(rz);
				}
			}
		  }

		  //Add the values to the descriptor vector
		  desc[count++] = dx;
		  desc[count++] = dy;
		  desc[count++] = dz;
		  desc[count++] = mdx;
		  desc[count++] = mdy;
		  desc[count++] = mdz;
		  len += (dx*dx + dy*dy + dz*dz + mdx*mdx + mdy*mdy + mdz*mdz);
		  k+=radius;

		}
      j += radius;
    }
    i += radius;
  }

  if ( !normalize ) return;

  //Convert to Unit Vector
  len = sqrt(len);
  for(int i = 0; i < ipt->size; ++i)
    desc[i] /= len;
}


//! Get the raw descriptor.
//! Modified deawDescriptor(int id)
void Surf::getRawDescriptor(int id, int radius)
{
  int y, x, z, sample_x, sample_y, sample_z, count=0;
  int i, j, k;
  double scale;
  double rx, ry, rz;

  Ipoint *ipt = &ipts[id];
  scale = ipt->scale;
  x = fRound(ipt->x);
  y = fRound(ipt->y);
  z = fRound(ipt->z);
  ipt->allocate( 3 * 8 * radius * radius * radius );
  float *desc = ipt->descriptor;

  i = -radius; //-10

  //Calculate descriptor for this interest point
  while(i < radius) // 10
  {
    j = -radius; //-10

    while(j < radius) //10
    {
		k = -radius;

		while(k < radius)
		{

		  for (int u = i; u < i + radius; ++u)
		  {
			for (int v = j; v < j + radius; ++v)
			{
				for (int w = k; w < k + radius; ++w)
				{
				  //Get coords of sample point on the rotated axis
				  sample_x = fRound(x + u*scale); // x + k
				  sample_y = fRound(y + v*scale); // y + l
				  sample_z = fRound(z + w*scale);

				  rx = haarXOptim(sample_x, sample_y, sample_z, 2*fRound(scale));
				  ry = haarYOptim(sample_x, sample_y, sample_z, 2*fRound(scale));
				  rz = haarZOptim(sample_x, sample_y, sample_z, 2*fRound(scale));

				  //Add the values to the descriptor vector
				  desc[count++] = rx;
				  desc[count++] = ry;
				  desc[count++] = rz;

				}
			}
		  }

		k+=radius;

		}
      j += radius;
    }
    i += radius;
  }

}

//-------------------------------------------------------

//! Calculate the value of the 2d gaussian at x,y,z (denormalized)
inline float Surf::gaussian(float x, float y, float z, float sig)
{
  return 1.0f/(sig*sig*sig)*exp( -(x*x+y*y+z*z)/(2.0f*sig*sig));
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in x direction
inline float Surf::haarX(int x, int y, int z, int s)
{
  return BoxIntegral(img, x, 		y-s/2,	z-s/2,	s/2, s, s)
    -1 * BoxIntegral(img, x-s/2, 	y-s/2,	z-s/2, 	s/2, s, s);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in y direction
inline float Surf::haarY(int x, int y, int z, int s)
{
  return BoxIntegral(img, x-s/2, 	y,		z-s/2,	s, s/2, s)
    -1 * BoxIntegral(img, x-s/2, 	y-s/2,	z-s/2, 	s, s/2, s);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in y direction
inline float Surf::haarZ(int x, int y, int z, int s)
{
  return BoxIntegral(img, x-s/2, 	y-s/2,	z,		s, s, s/2)
    - BoxIntegral(img, x-s/2, 	y-s/2,	z-s/2, 	s, s, s/2);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in x direction
inline float Surf::haarXOptim(int x, int y, int z, int s)
{
  
  return (long long)BoxIntegralOptim(img, x, y-s/2, z-s/2, s/2, s, s, incs, dims, pin_tegral)
  - (long long) BoxIntegralOptim(img, x-s/2, y-s/2, z-s/2, s/2, s, s, incs, dims, pin_tegral);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in y direction
inline float Surf::haarYOptim(int x, int y, int z, int s)
{
  return (long long)BoxIntegralOptim(img, x-s/2, 	y,		z-s/2,	s, s/2, s, incs, dims, pin_tegral)
    - (long long)BoxIntegralOptim(img, x-s/2, 	y-s/2,	z-s/2, 	s, s/2, s, incs, dims, pin_tegral);
}

//-------------------------------------------------------

//! Calculate Haar wavelet responses in y direction
inline float Surf::haarZOptim(int x, int y, int z, int s)
{
  return (long long)BoxIntegralOptim(img, x-s/2, 	y-s/2,	z,		s, s, s/2, incs, dims, pin_tegral)
    - (long long) BoxIntegralOptim(img, x-s/2, 	y-s/2,	z-s/2, 	s, s, s/2, incs, dims, pin_tegral);
}
