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

#ifndef FASTHESSIAN_H
#define FASTHESSIAN_H

#include <vector>

#include <vtkImageData.h>
#include <opencv2/opencv.hpp>
#include "ipoint.h"
#include <vtkMutexLock.h>

#define NONE_SCALE	0 //0b00
#define FIRST_SCALE 1 //0b01
#define LAST_SCALE  2 //0b10

class ResponseLayer;


//Maximum value
static const int OCTAVES = 5;
static const int INTERVALS = 4;
static const float THRES = 0.0004f;
static const int INIT_SAMPLE = 2;


class FastHessian {

  public:

    //! Constructor without image
    FastHessian(std::vector<Ipoint> &ipts,
				const int octaves = OCTAVES,
                const int intervals = INTERVALS,
                const int init_sample = INIT_SAMPLE,
                const float thres = THRES);

    //! Constructor with image
    FastHessian(vtkImageData *img,
				std::vector<Ipoint> &ipts,
                const int octaves = OCTAVES,
                const int intervals = INTERVALS,
                const int init_sample = INIT_SAMPLE,
                const float thres = THRES);

    //! Destructor
    ~FastHessian();

    //! Save the parameters
    void saveParameters(const int octaves,
                        const int intervals,
                        const int init_sample,
                        const float thres);

    //! Set or re-set the integral image source
    void setIntImage(vtkImageData *img);

    void setMask(vtkImageData *Mask);

    void WriteResponseMap();

    //! Find the image features and write into vector of features
    void getIpoints();

  private:

    //---------------- Private Functions -----------------//

    //! Build map of DoH responses
    void buildResponseMap();

    //! Calculate DoH responses for supplied layer
    void buildResponseLayer(int i, int decile);


    void EigenValue(float Cxx, float Cyy, float Czz, float Cxy, float Cyz, float Cxz, float c[3]);


    //! 3x3x3x3 Extrema test
    int isExtremum(int r, int c, int d, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, int param);
    int isCornerExtremum(int r, int c, int d, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, int param);

    //! Interpolation functions - adapted from Lowe's SIFT implementation
    void interpolateExtremum(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, bool corner);
    void interpolateStep(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                          double* xX, double* xY, double* xZ, double* xS, bool corner);
    CvMat* deriv4D(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, bool corner);

    CvMat* hessian4D(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, bool corner);

	static VTK_THREAD_RETURN_TYPE ThreadedRL (void *arg);

	vtkMutexLock* Mutex;
	int Current_thread_ID;
	int Current_decile;


	void FittingQuadric(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                                  double* xX, double* xY, double* xZ, double* xS , bool corner);

    //---------------- Private Variables -----------------//

    //! Pointer to the integral Image, and its attributes
    vtkImageData *img;
    int i_width, i_height, i_depth;

    //! Reference to vector of features passed from outside
    std::vector<Ipoint> &ipts;

	vtkImageData *Mask;


    //! Response stack of determinant of hessian values
    std::vector<ResponseLayer *> responseMap;

    //! Number of Octaves
    int octaves;

    //! Number of Intervals per octave
    int intervals;

    //! Initial sampling step for Ipoint detection
    int init_sample;

    //! Threshold value for blob resonses
    float thresh;

	inline void print(CvMat *input);
};


#endif
