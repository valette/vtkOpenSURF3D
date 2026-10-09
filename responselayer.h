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

#ifndef RESPONSELAYER_H
#define RESPONSELAYER_H

#include <cstring>
#include <cassert>
#include <cstddef>

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

  ResponseLayer(const ResponseLayer&) = delete;
  ResponseLayer& operator=(const ResponseLayer&) = delete;

  ResponseLayer(int width, int height, int depth, int step, int filter)
  {
    assert(width > 0 && height > 0 && depth > 0);

    this->width = width;
    this->height = height;
    this->depth = depth;
    this->step = step;
    this->filter = filter;

    size_t total = static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(depth);
    responses = new float[total]();
    cornerResponses = new float[total]();
    laplacian = new unsigned char[total]();
    isblob    = new bool[total]();
  }

  ~ResponseLayer()
  {
    delete [] responses;
    delete [] cornerResponses;
    delete [] laplacian;
    delete [] isblob;
  }

  inline unsigned char getLaplacian(unsigned int row, unsigned int column, unsigned int layer)
  {
    return laplacian[column + row * width + layer * width * height];
  }

  inline unsigned char getLaplacian(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;
    return laplacian[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

  inline float getResponse(unsigned int row, unsigned int column, unsigned int layer)
  {
    return responses[column + row * width + layer * width * height];
  }

  inline float getResponse(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;
    return responses[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

  inline float getCornerResponse(unsigned int row, unsigned int column, unsigned int layer)
  {
    return cornerResponses[column + row * width + layer * width * height];
  }

  inline float getCornerResponse(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;
    return cornerResponses[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

  inline bool getIsblob(unsigned int row, unsigned int column, unsigned int layer)
  {
    return isblob[column + row * width + layer * width * height];
  }

  inline bool getIsblob(unsigned int row, unsigned int column, unsigned int layer, ResponseLayer *src)
  {
    int scale = this->width / src->width;
    return isblob[(scale*column) + (scale*row) * width + (scale*layer) * width * height];
  }

};

#endif // RESPONSELAYER_H
