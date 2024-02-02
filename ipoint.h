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


//-------------------------------------------------------

class Ipoint; // Pre-declaration
typedef std::vector<Ipoint> IpVec;
typedef std::vector<std::pair<int, int> > IpPairVec;

//-------------------------------------------------------


class Ipoint {

public:

  //! Constructor
  Ipoint() : response(0), laplacian(0), scale(0) {};

  //! Gets the distance in descriptor space between Ipoints
  float operator-(const Ipoint &rhs)
  {
    float sum=0.f;
    for(int i=0; i < this->descriptor.size(); ++i)
      sum += (this->descriptor[i] - rhs.descriptor[i])*(this->descriptor[i] - rhs.descriptor[i]);
    return sqrt(sum);
  };

  void allocate( int size ) {
    this->descriptor.resize( size );
  }

  //! Coordinates of the detected interest point
  float x, y, z;

  //! Detected scale
  float scale;
  
  //! Response
  float response;

  //! Sign of laplacian for fast matching purposes
  int laplacian;

  //! Vector of descriptor components
  std::vector< float > descriptor;

  //! Point id
//  int id;
};

//-------------------------------------------------------


#endif
