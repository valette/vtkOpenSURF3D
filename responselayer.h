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

#include <memory.h>
#include <assert.h>

// x = width =  column
// y = height = row
// z = depth =  layer
// x + y * width + z * width * height
// column + row * width + layer * width * height

class ResponseLayer
{
public:

  int width, height, depth, step, filter;
  float *responses;
  unsigned char *laplacian;
  float *cornerResponses;
    // 0 : dark blob,
    // 1 ; light blob,
    // 2 : corner
  bool *isblob;

  ResponseLayer(int width, int height, int depth, int step, int filter)
  {
    assert(width > 0 && height > 0 && depth > 0);

    this->width = width;
    this->height = height;
    this->depth = depth;
    this->step = step;
    this->filter = filter;

    responses = new float[width*height*depth];
    cornerResponses = new float[width*height*depth];
    laplacian = new unsigned char[width*height*depth];
    isblob    = new bool[width*height*depth];

    memset(responses,0,sizeof(float)*width*height*depth);
    //memset(laplacian,0,sizeof(unsigned char)*width*height*depth);
    //memset(isblob,0,sizeof(bool)*width*height*depth);
  }

  ~ResponseLayer()
  {
    if (responses) delete [] responses;
    if (cornerResponses) delete [] cornerResponses;
    if (laplacian) delete [] laplacian;
    if (isblob) delete [] isblob;
  }

  inline unsigned char getLaplacian(unsigned int row, unsigned int column, unsigned int layer)
  {
	//assert(row < height && width > column && depth >  layer);
    return laplacian[column + row * width + layer * width * height];
  }

  inline unsigned char getLaplacian(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;

    //assert(scale*row < height && width > scale*column && depth >  scale*layer);

    return laplacian[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

  inline float getResponse(unsigned int row, unsigned int column, unsigned int layer)
  {
	  	//assert(row < height && width > column && depth >  layer);

    return responses[column + row * width + layer * width * height];
  }

  inline float getResponse(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;

    //assert(scale*row < height && width > scale*column && depth >  scale*layer);

    return responses[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

  inline float getCornerResponse(unsigned int row, unsigned int column, unsigned int layer)
  {
      //assert(row < height && width > column && depth >  layer);

    return cornerResponses[column + row * width + layer * width * height];
  }

  inline float getCornerResponse(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;

    //assert(scale*row < height && width > scale*column && depth >  scale*layer);

    return cornerResponses[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

  inline float getIsblob(unsigned int row, unsigned int column, unsigned int layer)
  {
	//assert(row < height && width > column && depth >  layer);

    return isblob[column + row * width + layer * width * height];
  }

  inline float getIsblob(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;

    //assert(scale*row < height && width > scale*column && depth >  scale*layer);

    return isblob[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

};
