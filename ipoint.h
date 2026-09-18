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
#include <cmath>
#include <algorithm>

//-------------------------------------------------------

class Ipoint; // Pre-declaration
typedef std::vector<Ipoint> IpVec;
typedef std::vector<std::pair<int, int> > IpPairVec;

//-------------------------------------------------------

class Ipoint {

public:

  //! Constructor
  Ipoint() : x(0.0f), y(0.0f), z(0.0f), scale(0.0f), response(0.0f), laplacian(0) {}

  //! Gets the distance in descriptor space between Ipoints
  float operator-(const Ipoint &rhs) const
  {
    float sum = 0.0f;
    size_t minSize = std::min(this->descriptor.size(), rhs.descriptor.size());
    for(size_t i = 0; i < minSize; ++i)
    {
      float diff = this->descriptor[i] - rhs.descriptor[i];
      sum += diff * diff;
    }
    return std::sqrt(sum);
  }

  void allocate( int size ) {
    this->descriptor.resize( size );
  }

  //! Coordinates of the detected interest point
  float x;
  float y;
  float z;

  //! Detected scale
  float scale;
  
  //! Response
  float response;

  //! Sign of laplacian for fast matching purposes
  int laplacian;

  //! Vector of descriptor components
  std::vector< float > descriptor;

};

//-------------------------------------------------------


#endif
