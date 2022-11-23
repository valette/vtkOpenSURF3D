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

#ifndef IPOINT_H
#define IPOINT_H

#include <vector>
#include <math.h>
#include <limits>


//-------------------------------------------------------

class Ipoint; // Pre-declaration
typedef std::vector<Ipoint> IpVec;
typedef std::vector<std::pair<Ipoint, Ipoint> > IpPairVec;

//-------------------------------------------------------

//! Ipoint operations
void getMatches(IpVec &ipts1, IpVec &ipts2, IpPairVec &matches);
//-------------------------------------------------------

class Ipoint {

public:

  //! Destructor
  ~Ipoint() {
 //   if (descriptor) delete descriptor;
  };

  //! Constructor
  Ipoint() : orientation(0), descriptor(0), response(0), laplacian(0), scale(0), size(0) {};

  //! Gets the distance in descriptor space between Ipoints
  float operator-(const Ipoint &rhs)
  {
    float sum=0.f;
    for(int i=0; i < this->size; ++i)
      sum += (this->descriptor[i] - rhs.descriptor[i])*(this->descriptor[i] - rhs.descriptor[i]);
    return sqrt(sum);
  };

  void allocate( int size ) {
	  this->size = size;
	  this->descriptor = new float[ size ];
  }

  //! Coordinates of the detected interest point
  float x, y, z;

  //! Detected scale
  float scale;
  
  //! Response
  float response;

  //! Orientation measured anti-clockwise from +ve x-axis
  float orientation;

  //! Sign of laplacian for fast matching purposes
  int laplacian;

  //! Size of descriptor
  int size;

  //! Vector of descriptor components
  float *descriptor;

  //! Point id
  int id;
};

//-------------------------------------------------------


#endif
