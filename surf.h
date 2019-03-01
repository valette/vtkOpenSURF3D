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

#ifndef SURF_H
#define SURF_H

#include <opencv2/opencv.hpp>
#include "ipoint.h"
#include "integral.h"

#include <vector>

typedef unsigned long long int ulli;


class Surf {
  
  public:
    
    //! Standard Constructor (img is an integral image)
    Surf(vtkImageData *img, std::vector<Ipoint> &ipts);

    //! Describe all features in the supplied vector
    void getDescriptors( int radius = 5 );

    //! Describe all raw features in the supplied vector
    void getRawDescriptors( int radius = 5 );
  
  private:
    
    //---------------- Private Functions -----------------//
  
    static VTK_THREAD_RETURN_TYPE ThreadedDesc (void *arg);
    
    //! Assign the current Ipoint an orientation
    //void getOrientation();
    
    //! Get the descriptor. See Agrawal ECCV 08
    void getDescriptor(int id, int radius);

    //! Get the raw descriptor.
    void getRawDescriptor(int id, int radius);

    //! Calculate the value of the 2d gaussian at x,y
    inline float gaussian(float x, float y, float z, float sig);
    //inline float gaussian(float x, float y, float sig);

    //! Calculate Haar wavelet responses in x and y directions
    inline float haarX(int x, int y, int z, int s);
    inline float haarY(int x, int y, int z, int s);
    inline float haarZ(int x, int y, int z, int s);

	inline float haarXOptim(int x, int y, int z, int s);
	inline float haarYOptim(int x, int y, int z, int s);
	inline float haarZOptim(int x, int y, int z, int s);


    //! Get the angle from the +ve x-axis of the vector given by [X Y]
    //float getAngle(float X, float Y);


    //---------------- Private Variables -----------------//

    //! Integral image where Ipoints have been detected
    vtkImageData *img;

    //! Ipoints vector
    IpVec &ipts;
    
    int* dims;
	long long int* incs;
	unsigned long long int* pin_tegral;

};


//! Round float to nearest integer
inline int fRound(float flt)
{
  return (int) floor(flt+0.5f);
}


#endif
