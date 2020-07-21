#include "integral.h"
#include "ipoint.h"
//#include "utils.h"

#include <opencv2/opencv.hpp>

#include <vector>

#include "responselayer.h"
#include "fasthessian.h"

#include <vtkImageData.h>
#include <vtkMultiThreader.h>
#include <vtkTimerLog.h>
#include <vtkMetaImageWriter.h>
#include <vtkType.h>

#include <math.h>


using namespace std;

int nb_pts = 0;
int nb_corner_pts = 0;

//-------------------------------------------------------

//! Constructor without image
FastHessian::FastHessian(std::vector<Ipoint> &ipts,
						 const int octaves, const int intervals, const int init_sample,
                         const float thresh)
                         : ipts(ipts), i_width(0), i_height(0), i_depth(0), Current_thread_ID(0), Current_decile(0), Mask(0)
{
  // Save parameter set
  saveParameters(octaves, intervals, init_sample, thresh);
}

//-------------------------------------------------------

//! Constructor with image
FastHessian::FastHessian(vtkImageData *img,
						 std::vector<Ipoint> &ipts,
                         const int octaves, const int intervals, const int init_sample,
                         const float thresh)
                         : ipts(ipts), i_width(0), i_height(0), i_depth(0), Current_thread_ID(0), Current_decile(0), Mask(0)
{
  // Save parameter set
  saveParameters(octaves, intervals, init_sample, thresh);

  // Set the current image
  setIntImage(img);

  int* dims = img->GetDimensions();
  int minDim = std::min({dims[0], dims[1], dims[2]});

  if (minDim < 51*2)
  {
    this->octaves = min(1, this->octaves);
  }
  else if (minDim < 99*2)
  {
    this->octaves = min(2, this->octaves);
  }
  else if (minDim < 195*2)
  {
    this->octaves = min(3, this->octaves);
  }
  else if (minDim < 387*2)
  {
    this->octaves = min(4, this->octaves);
  }
}

//-------------------------------------------------------

FastHessian::~FastHessian()
{
  for (unsigned int i = 0; i < responseMap.size(); ++i)
  {
    delete responseMap[i];
  }
}

//-------------------------------------------------------

//! Save the parameters
void FastHessian::saveParameters(const int octaves, const int intervals,
                                 const int init_sample, const float thresh)
{
  // Initialise variables with bounds-checked values
  this->octaves =
    (octaves > 0 && octaves <= 4 ? octaves : OCTAVES);
  this->intervals =
    (intervals > 0 && intervals <= 4 ? intervals : INTERVALS);
  this->init_sample =
    (init_sample > 0 && init_sample <= 6 ? init_sample : INIT_SAMPLE);
  this->thresh = (thresh >= 0 ? thresh : THRES);
}


//-------------------------------------------------------

//! Set or re-set the integral image source
void FastHessian::setIntImage(vtkImageData *img)
{
  // Change the source image
  this->img = img;

  int* dims = img->GetDimensions();

  i_height = dims[1];
  i_width = dims[0];
  i_depth = dims[2];
}

//-------------------------------------------------------

//! Set or re-set the integral image source
void FastHessian::setMask(vtkImageData *Mask)
{
  // Change the mask image
  this->Mask = Mask;
}


//-------------------------------------------------------

//! Find the image features and write into vector of features
void FastHessian::getIpoints()
{
  // filter index map
  static const int filter_map [OCTAVES][INTERVALS] = {{0,1,2,3}, {1,3,4,5}, {3,5,6,7}, {5,7,8,9}, {7,9,10,11}};

  // Clear the vector of exisiting ipts
  //ipts.clear();


	vtkTimerLog *Timer = vtkTimerLog::New();

	Timer->StartTimer();
  // Build the response map
  buildResponseMap();

	Timer->StopTimer();
	cout << "response Map built (" << Timer->GetElapsedTime() << "s)" << endl;


	Timer->StartTimer();

	int layer_min = filter_map[0][0];
	int layer_max = filter_map[octaves][2];

  // Get the response layers
  ResponseLayer *b, *m, *t;
  for (int o = 0; o < octaves; ++o) for (int i = 0; i <= 1; ++i)
  {

			b = responseMap.at(filter_map[o][i]);
			m = responseMap.at(filter_map[o][i+1]);
			t = responseMap.at(filter_map[o][i+2]);

			int param;
			if (o == 0 && i == 0)
				param = FIRST_SCALE;
			else if (o == (octaves-1) && i == 1)
				param = LAST_SCALE;
			else
				param = NONE_SCALE;

	//int limit = std::max((float)ceil((float)((m->filter - 1) / 2)/(float)m->step)+2, ((float)(0.6666*m->filter + 2.5 * m->step)/m->step+2));
	int limit_m = ceil((ceil((float)(m->filter + 1)/(float)m->step/2.0)+1)* (float)m->step/(float)t->step);
	int limit_b = ceil((ceil((float)(b->filter + 1)/(float)b->step/2.0)+1)* (float)b->step/(float)t->step);
	int limit_t = ceil((float)(t->filter + 1)/(float)t->step/2.0)+1;

	int limit = max(limit_m, limit_b);
	limit = max(limit, limit_t);

	limit = max((float)limit, (float)(0.6666*t->filter + 2.5 * m->step)/(float)m->step + 2);


    // loop over middle response layer at density of the most
    // sparse layer (always top), to find maxima across scale and space
    for (int r = limit; r < t->height-limit; ++r)
      for (int c = limit; c < t->width-limit; ++c)
		  for (int d = limit; d < t->depth-limit; ++d)
		  {

			int LimSupScale  = std::max((float)ceil((float)((t->filter - 1) / 2 +2 )/(float)t->step), (float)(0.6666*t->filter + 2.5 * t->step +2)/t->step );
			int LimDownScale = std::max((float)ceil((float)((b->filter - 1) / 2 +2 )/(float)b->step), (float)(0.6666*b->filter + 2.5 * b->step  +2)/b->step);

			if ( r < LimSupScale || r >= t->height- LimSupScale ||
				 c < LimSupScale || c >= t->width - LimSupScale ||
				 d < LimSupScale || d >= t->depth - LimSupScale )
					param &= LAST_SCALE;

			if ( r*(int)(t->width/b->width) < LimDownScale || r*(int)(t->width/b->width) >= b->height- LimDownScale ||
				 c*(int)(t->width/b->width) < LimDownScale || c*(int)(t->width/b->width) >= b->width - LimDownScale ||
				 d*(int)(t->width/b->width) < LimDownScale || d*(int)(t->width/b->width) >= b->depth - LimDownScale )
					param &= FIRST_SCALE;

			//param = 0;
			//cout << "isExtremum(" << r << "," << c << "," << d << "," << o << "," << i << "," << param  << ");" << endl;

			if (isExtremum(r, c, d, t, m, b, param))
			{
			  interpolateExtremum(d, r, c, t, m, b, false);
			}
#ifdef EXTRACT_CORNER
			if (isCornerExtremum(r, c, d, t, m, b, param))
			{
			  interpolateExtremum(d, r, c, t, m, b, true);
			}
#endif

		}
	}

	Timer->StopTimer();
	cout << "find maxima across scale and space (" << Timer->GetElapsedTime() << "s)" <<endl;


	//Remove OutMask Point
	if (Mask != 0 && false)
	{
		cout << ipts.size() << " Ipoints before removal" << endl;

		vector<Ipoint>::iterator it = ipts.begin();

		while (it!=ipts.end())
		{
			bool flag = false;

			for (float i = -1 ; i < 2 ; i++)
			{
				for (float j = -1 ; j < 2 ; j++)
				{
					for (float k = -1 ; k < 2 ; k++)
					{
						if ( * static_cast<unsigned char*>(Mask->GetScalarPointer(
								it->x + i * 0 * it->scale,
								it->y + j * 0 * it->scale,
								it->z + k * 0 * it->scale)) == 0)
								{
							it = ipts.erase(it);
							flag = true;
							break;
						}
					}
					if (flag) break;
				}
				if (flag) break;
			}
			if (!flag) it++;
		}
			cout << ipts.size() << " Ipoints after removal" << endl;
	}
	else {
		cout << " Ipoints : " << ipts.size() << endl;
		cout << nb_corner_pts << "\tcorner" << endl;
		cout << nb_pts << "\tblobs" << endl;
	}



}

//-------------------------------------------------------

//! Build map of DoH responses
void FastHessian::buildResponseMap()
{
  // Calculate responses for the first 4 octaves:
  // Oct1: 9,  15, 21, 27
  // Oct2: 15, 27, 39, 51
  // Oct3: 27, 51, 75, 99
  // Oct4: 51, 99, 147,195
  // Oct5: 99, 195,291,387

  // Deallocate memory and clear any existing response layers
  for(unsigned int i = 0; i < responseMap.size(); ++i)
    delete responseMap[i];
  responseMap.clear();

  // Get image attributes
  int w = (i_width / init_sample);
  int h = (i_height / init_sample);
  int d = (i_depth / init_sample);
  int s = (init_sample);

  // Calculate approximated determinant of hessian values
  if (octaves >= 1)
  {
    responseMap.push_back(new ResponseLayer(w,   h,   d,   s,   9));
    responseMap.push_back(new ResponseLayer(w,   h,   d,   s,   15));
    responseMap.push_back(new ResponseLayer(w,   h,   d,   s,   21));
    responseMap.push_back(new ResponseLayer(w,   h,   d,   s,   27));
  }

  if (octaves >= 2)
  {
    responseMap.push_back(new ResponseLayer(w/2, h/2, d/2, s*2, 39));
    responseMap.push_back(new ResponseLayer(w/2, h/2, d/2, s*2, 51));
  }

  if (octaves >= 3)
  {
    responseMap.push_back(new ResponseLayer(w/4, h/4, d/4, s*4, 75));
    responseMap.push_back(new ResponseLayer(w/4, h/4, d/4, s*4, 99));
  }

  if (octaves >= 4)
  {
    responseMap.push_back(new ResponseLayer(w/8, h/8, d/8, s*8, 147));
    responseMap.push_back(new ResponseLayer(w/8, h/8, d/8, s*8, 195));
  }

  if (octaves >= 5)
  {
    responseMap.push_back(new ResponseLayer(w/16, h/16, d/16, s*16, 291));
    responseMap.push_back(new ResponseLayer(w/16, h/16, d/16, s*16, 387));
  }

  // Extract responses from the image

  for (unsigned int i = 0; i < responseMap.size(); ++i)
	buildResponseLayer(i, 0);

}

//! Calculate DoH responses for supplied layer
void FastHessian::buildResponseLayer(int id, int decile)
{

  ResponseLayer* rl = responseMap[id];

  float *responses = rl->responses;         // response storage
  unsigned char *laplacian = rl->laplacian; // laplacian sign storage
  bool *isblob = rl->isblob;					// isblob storage
	float *cornerResponses = rl->cornerResponses;
  int step = rl->step;                      // step size for this filter
  int b = (rl->filter - 1) / 2;             // border for this filter
  int l = rl->filter / 3;                   // lobe for this filter (filter size / 3)
  int w = rl->filter;                       // filter size
  float inverse_volume_cube = 1.f/(((float)w*(float)w*(float)w)*((float)w*(float)w*(float)w)*((float)w*(float)w*(float)w));           // normalisation factor

  //Accesseur rapide aux voxels
  int* dim = img->GetDimensions();
  long long int* incs = img->GetIncrements();
  unsigned long long int* pin_tegral = static_cast<ulli*>(img->GetScalarPointer());

  //Limits pour etre toujours à l'interieur de l'image (Box integral & Descripteurs)
  // 5*(0.1333 * m->filter + xS * filterStep) (xS < 0.5)

  int limit = ceil((float)(rl->filter + 1)/(float)step/2)+1; //std::max((float)(0.6666*w + 2.5 * step+1)/step,

  //Limit calculé dans l'espace du response layer (en terme d'indice)

  //if (decile+2*limit >= rl->width)
  //	return;

  float tmax = 5;
  float curv = (2 * tmax + 1.0)*(2 * tmax + 1.0)*(2 * tmax + 1.0) / (tmax * tmax);



  //for(int x, y, z, index, ax = decile+limit; ax < rl->width-limit; ax += GRANULARITE
 #pragma omp parallel for shared(rl, isblob, responses, laplacian)
  for(int az = limit; az < rl->depth-limit; ++az)
  {
    for(int ay = limit; ay < rl->height-limit; ++ay)
    {
  		for(int ax = limit; ax < rl->width-limit; ++ax)
		{
			//cout << "response " << ax << " " << ay << " " << az <<  ":" << id << endl;

			int index = ax + ay * rl->width + az * rl->width * rl->height;
			  // get the image coordinates
			  int x = ax * step;
			  int y = ay * step;
			  int z = az * step;
			  /*
			  if (x+b >= dim[0] || y+b >= dim[1] || z+b >= dim[2])
			  {
				cout << "GAAA" << (bool)(x+b >= dim[0]) << (bool)(y+b >= dim[1]) << (bool)(z+b >= dim[2]) <<  endl;
				cout << "limit="<< limit << endl;
				cout << "dim[1]="<< dim[1] << endl;
				cout << "y="<< y << endl;
				cout << "b="<< b << endl;
				cout << "rl->height=" << rl->height << endl;
				cout << "step=" << step << endl;
			}
			  */

				float Dxx = (float) BoxIntegralOptim(img, 	x-b,	y-l+1,	z-l+1, 	w,		2*l-1,	2*l-1, 	incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x-l/2, 	y-l+1, 	z-l+1, 	l, 		2*l-1,	2*l-1,	incs, dim, pin_tegral)*3;

				float Dyy = (float) BoxIntegralOptim(img, 	x-l+1,	y-b,	z-l+1, 	2*l-1,	w,		2*l-1, 	incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x-l+1, 	y-l/2, 	z-l+1, 	2*l-1, 	l,		2*l-1, 	incs, dim, pin_tegral)*3;

				float Dzz = (float) BoxIntegralOptim(img, 	x-l+1,	y-l+1,	z-b, 	2*l-1,	2*l-1,	w, 		incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x-l+1, 	y-l+1, 	z-l/2, 	2*l-1, 	2*l-1,	l, 		incs, dim, pin_tegral)*3;

			  float Dxy = + (float) BoxIntegralOptim(img, 	x - l, 	y - l,	z-l+1,	l, 		l, 		2*l-1, 	incs, dim, pin_tegral)
					+ (float) BoxIntegralOptim(img, 	x + 1, 	y + 1, 	z-l+1,	l, 		l, 		2*l-1, 	incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x - l, 	y + 1, 	z-l+1,	l, 		l, 		2*l-1, 	incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x + 1, 	y - l, 	z-l+1,	l, 		l, 		2*l-1, 	incs, dim, pin_tegral);

			  float Dyz = + (float) BoxIntegralOptim(img, 	x-l+1, 	y - l,	z - l,	2*l-1, 	l, 		l, 		incs, dim, pin_tegral)
					+ (float) BoxIntegralOptim(img, 	x-l+1, 	y + 1,	z + 1, 	2*l-1, 	l, 		l, 		incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x-l+1, 	y - l,	z + 1, 	2*l-1, 	l, 		l, 		incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x-l+1, 	y + 1,	z - l, 	2*l-1, 	l, 		l, 		incs, dim, pin_tegral);

			  float Dxz = + (float) BoxIntegralOptim(img, 	x - l, 	y-l+1,	z - l,	l, 		2*l-1,	l, 		incs, dim, pin_tegral)
					+ (float) BoxIntegralOptim(img, 	x + 1, 	y-l+1, 	z + 1,	l, 		2*l-1,	l, 		incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x - l, 	y-l+1, 	z + 1,	l, 		2*l-1,	l, 		incs, dim, pin_tegral)
					- (float) BoxIntegralOptim(img, 	x + 1, 	y-l+1, 	z - l,	l, 		2*l-1,	l, 		incs, dim, pin_tegral);
#ifdef EXTRACT_CORNER

				float Ix = + (float) BoxIntegralOptim(img, 	x - l, 	y-l+1,	z-l+1,	l, 		2*l-1,	2*l-1, 		incs, dim, pin_tegral)
					+ (float) BoxIntegralOptim(img, 	x + 1, 	y-l+1, 	z-l+1,	l, 		2*l-1,	2*l-1, 		incs, dim, pin_tegral);

				float Iy = + (float) BoxIntegralOptim(img, 	x-l+1, 	y - l,	z-l+1,	l*l-1, 		l,	2*l-1, 		incs, dim, pin_tegral)
					+ (float) BoxIntegralOptim(img, 	x-l+1, 	y + 1, 	z-l+1,	l*l-1, 		l,	2*l-1, 		incs, dim, pin_tegral);

				float Iz = + (float) BoxIntegralOptim(img, 	x-l+1, 	y-l+1,	z - l,	l*l-1, 		2*l-1,	l, 		incs, dim, pin_tegral)
					+ (float) BoxIntegralOptim(img, 	x-l+1, 	y-l+1, 	z + 1,	l*l-1, 		2*l-1,	l, 		incs, dim, pin_tegral);
#endif
		//Quand rien dans la dim : dim = r-l+1 		size = 2*l-1

				float Sdet2p = Dyy*Dzz + Dxx*Dyy + Dxx*Dzz - 0.8330f*(Dxy*Dxy + Dxz*Dxz + Dyz*Dyz);
				float Trace  = Dxx + Dyy + Dzz;
				float Det	 = Dxx*Dyy*Dzz + 2.0*Dxy*Dyz*Dxz*0.7603f - Dxx*Dyz*Dyz*0.8330f - Dyy*Dxz*Dxz*0.8330f - Dzz*Dxy*Dxy*0.8330f;


				isblob[index] =  (Sdet2p > 0) && (Trace*Det>0); //&& (Trace*Trace*Trace/Det < curv);

			  // Get the determinant normalized of hessian response & laplacian sign

			  responses[index] = abs(Det*inverse_volume_cube);

			  laplacian[index] = (Dxx + Dyy + Dzz >= 0 ? 1 : 0);

#ifdef EXTRACT_CORNER
/*
				float CornerDet =+ Ix*Ix * (Iy*Iy * Iz*Iz - Iz*Iy * Iy*Iz)
                         - Ix*Iy * (Iy*Ix * Iz*Iz - Iy*Iz * Iz*Ix)
                         + Ix*Iz * (Iy*Ix * Iz*Iy - Iy*Iy * Iz*Ix);
*/
				float eigs[3];
				float a[] = {
				Ix*Ix , Ix*Iy , Ix*Iz ,
				Ix*Iy , Iy*Iy  ,Iy*Iz  ,
				Ix*Iz , Iy*Iz ,  Iz*Iz };

				CvMat mat = cvMat(3,3,CV_32FC1, a);
				CvMat* evec  = cvCreateMat(3,3,CV_32FC1);
				CvMat* eval  = cvCreateMat(3,1,CV_32FC1);

				cvZero(evec);
				cvZero(eval);

				cvEigenVV(&mat, evec, eval, DBL_EPSILON, -1, -1);

			  //EigenValue(Ix*Ix, Iy*Iy, Iz*Iz, Ix*Iy, Iy*Iz, Ix*Iz, eigs);
				//cornerResponses[index] = std::min(std::min(abs(eigs[0]), abs(eigs[1])), abs(eigs[2])) * inverse_volume_cube;
				cornerResponses[index] = std::min(std::min(cvGet2D(&mat,0,0).val[0], cvGet2D(&mat,0,1).val[0]), cvGet2D(&mat,0,2).val[0]) * 1.f/((float)w*(float)w*(float)w*(float)w*(float)w*(float)w);
#endif
			}
		}
	}
}

//-------------------------------------------------------
// http://www.rohitab.com/discuss/topic/36251-c-svd-of-3x3-matrix/
void FastHessian::EigenValue(float Cxx, float Cyy, float Czz, float Cxy, float Cyz, float Cxz, float c[3]) {

	c[2] = -Cxx - Cyy - Czz;
	c[1] = (Cxx*Cyy + Czz*Cxx + Czz*Cyy -
					Cyz*Cyz - Cxz*Cxz - Cxy*Cxy)*0.5;
	c[0] = Cyz*Cyz*Cxx + Cxz*Cxz*Cyy + Cxy*Cxy*Czz -
					Cxx*Cyy*Czz - Cxy*Cyz*Cxz - Cxz*Cxy*Cyz;


	const double sq3d2 = 0.86602540378443864676, c2d3 = c[2]/3,
          c2sq = c[2]*c[2], Q = (3*c[1]-c2sq)/9,
          R = (c[2]*(9*c[1]-2*c2sq)-27*c[0])/54;
  double tmp, t, sint, cost;

  if (Q < 0) {
    tmp = 2*sqrt(-Q);
    t = acos(R/sqrt(-Q*Q*Q))/3;
    cost = tmp*cos(t);
    sint = tmp*sin(t);

    c[0] = cost - c2d3;

    cost = -0.5*cost - c2d3;
    sint = sq3d2*sint;

    c[1] = cost - sint;
    c[2] = cost + sint;
  }
  else {
    tmp = cbrt(R);
    c[0] = -c2d3 + 2*tmp;
    c[1] = c[2] = -c2d3 - tmp;
  }
}

//! Non Maximal Suppression function
int FastHessian::isExtremum(int r, int c, int d, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, int param)
{

  // check the candidate point in the middle layer is above thresh
  float candidate = m->getResponse(r, c, d, t);

  if ( candidate < thresh || !m->getIsblob(r, c, d, t) )
  {
    return 0;
   }

  for (int rr = -1; rr <=1; ++rr)
    for (int cc = -1; cc <=1; ++cc)
		for (int dd = -1; dd <=1; ++dd)
		{
		  // if any response in 3x3x3x3 is greater candidate not maximum (only if isblob == true)
		  if (
			(t->getResponse(r+rr, c+cc, d+dd)    >= candidate && param != LAST_SCALE  && t->getIsblob(r+rr, c+cc, d+dd)) ||
			(m->getResponse(r+rr, c+cc, d+dd, t) >= candidate && (rr != 0 || cc != 0) && m->getIsblob(r+rr, c+cc, d+dd, t)) ||
			(b->getResponse(r+rr, c+cc, d+dd, t) >= candidate && param != FIRST_SCALE && b->getIsblob(r+rr, c+cc, d+dd, t))
			)
				return 0;
		}

  return true;
}

//! Non Maximal Suppression function
int FastHessian::isCornerExtremum(int r, int c, int d, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, int param)
{
  // check the candidate point in the middle layer is above thresh
  float candidate = m->getCornerResponse(r, c, d, t);

	//cout << "test extremum : " << candidate << endl;

  if ( candidate < thresh)
  {
    return 0;
   }

  for (int rr = -1; rr <=1; ++rr)
    for (int cc = -1; cc <=1; ++cc)
		for (int dd = -1; dd <=1; ++dd)
		{
		  // if any response in 3x3x3x3 is greater candidate not maximum (only if isblob == true)
		  if (
			(t->getCornerResponse(r+rr, c+cc, d+dd)    >= candidate && param != LAST_SCALE) ||
			(m->getCornerResponse(r+rr, c+cc, d+dd, t) >= candidate && (rr != 0 || cc != 0) ) ||
			(b->getCornerResponse(r+rr, c+cc, d+dd, t) >= candidate && param != FIRST_SCALE )
			)

				//cout << t->getCornerResponse(r+rr, c+cc, d+dd) << " " << m->getCornerResponse(r+rr, c+cc, d+dd, t) << " " << b->getCornerResponse(r+rr, c+cc, d+dd, t) << endl;
				return 0;
		}
		//cout << "yes ! " << endl;
  return true;
}

//-------------------------------------------------------
 int fsdi = 0;

//! Interpolate scale-space extrema to subpixel accuracy to form an image feature.
void FastHessian::interpolateExtremum(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, bool corner)
{
  // get the step distance between filters
  // check the middle filter is mid way between top and bottom
  int filterStep = (m->filter - b->filter);
  //assert(filterStep > 0 && t->filter - m->filter == m->filter - b->filter);

  // Get the offsets to the actual location of the extremum
  double xX = 0, xY = 0, xZ = 0, xS = 0;
  interpolateStep(d, r, c, t, m, b, &xX, &xY, &xZ, &xS, corner );

  //FittingQuadric( d, r, c, t, m, b,  &xX, &xY, &xZ, &xS );

  // If point is sufficiently close to the actual extremum
  if( fabs( xX ) < 1.0f  &&  fabs( xY ) < 1.0f  &&  fabs( xZ ) < 1.0f  &&  fabs( xS ) < 1.0f )
  {
    Ipoint ipt;
    ipt.x = static_cast<float>((c + xX) * t->step);
    ipt.y = static_cast<float>((r + xY) * t->step);
    ipt.z = static_cast<float>((d + xZ) * t->step);
    ipt.scale = static_cast<float>((0.1333f) * (m->filter + xS * filterStep));
#ifdef EXTRACT_CORNER
		if (corner) {
			nb_corner_pts++;
			ipt.laplacian = 2;
    	ipt.response = static_cast<float>(m->getCornerResponse(r,c,d,t));
		} else {
#endif
			nb_pts++;
			ipt.laplacian = static_cast<int>(m->getLaplacian(r,c,d,t));
			ipt.response = static_cast<float>(m->getResponse(r,c,d,t));
#ifdef EXTRACT_CORNER
		}
#endif
    ipts.push_back(ipt);
   }
}

//-------------------------------------------------------
//Attention différent de surf classique : les matrices sont dans l'ordre : x y z s
//! Performs one step of extremum interpolation.
void FastHessian::interpolateStep(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                                  double* xX, double* xY, double* xZ, double* xS, bool corner)
{
  double x[4] = { 0 };

  cv::Matx41f* dD = deriv4D( d, r, c, t, m, b, corner );
  cv::Matx44f* H = hessian4D( d, r, c, t, m, b, corner );
  cv::Matx44f H_inv, W_inv, V, U_t, tmp;
  cv::Matx41f W;

  W_inv = cv::Matx44f::zeros();

  cv::SVD::compute(*H, W, U_t, V);

	bool flag = false;

	W_inv( 0, 0 ) = 1 / W(0);

	for (int i = 1 ; i < 4 ; i++)
		if ( W(i) / W(0) < 0.001)
		{
			W_inv( i, i ) = 0;
		}
		else
		{
			W_inv( i, i ) = 1.0f / W( i);
			flag = true;
		}

//	cv::MatOp::matmul(W_inv, U_t, tmp);
	tmp = W_inv * U_t;
//	cv::MatOp::matmul(V, tmp, H_inv);
	H_inv = V *tmp;

  cv::Matx41f X = cv::Matx41f::zeros();


  //cvInvert( H, H_inv, CV_SVD );


//  cvInitMatHeader( &X, 4, 1, CV_64FC1, x, CV_AUTOSTEP );
//  cvGEMM( H_inv, dD, -1, NULL, 0, &X, 0 );

//  X = - H_inv.t()*dD;
  X = - H_inv* (*dD);


  delete dD;
  delete H;

  *xX = X(0);
  *xY = X(1);
  *xZ = X(2);
  *xS = X(3);
  /*
  if( fabs( *xX ) > 1.0f  ||  fabs( *xY ) > 1.0f  ||  fabs( *xZ ) > 1.0f  ||  fabs( *xS ) > 1.0f )
	{
		cout << "Points rejeté valeurs propres : " << log(cvmGet(W, 3, 3)/cvmGet(W, 0, 0)) << " flag : " << flag << endl;
	}
  */
}

//-------------------------------------------------------

//! Computes the partial derivatives in x, y, and scale of a pixel.
cv::Matx41f* FastHessian::deriv4D(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, bool corner)
{
  double dx, dy, ds, dz;

	if (corner) {
		dx = (m->getCornerResponse(r, c + 1, d, t) - m->getCornerResponse(r, c - 1,d,  t)) / 2.0;
	  dy = (m->getCornerResponse(r + 1, c, d, t) - m->getCornerResponse(r - 1, c,d,  t)) / 2.0;
	  dz = (m->getCornerResponse(r, c, d + 1, t) - m->getCornerResponse(r, c, d - 1,t)) / 2.0;
	  ds = (t->getCornerResponse(r, c, d) - 		 b->getCornerResponse(r, c, d, t)) / 2.0;
	} else {
		dx = (m->getResponse(r, c + 1, d, t) - m->getResponse(r, c - 1,d,  t)) / 2.0;
	  dy = (m->getResponse(r + 1, c, d, t) - m->getResponse(r - 1, c,d,  t)) / 2.0;
	  dz = (m->getResponse(r, c, d + 1, t) - m->getResponse(r, c, d - 1,t)) / 2.0;
	  ds = (t->getResponse(r, c, d) - 		 b->getResponse(r, c, d, t)) / 2.0;
	}

  cv::Matx41f* dI = new cv::Matx41f();
  (*dI)(0) = dx;
  (*dI)(1) = dy;
  (*dI)(2) = dz;
  (*dI)(3) = ds;

//  cvmSet( dI, 0, 0, dx );
//  cvmSet( dI, 1, 0, dy );
//  cvmSet( dI, 2, 0, dz );
//  cvmSet( dI, 3, 0, ds );

  return dI;
}


#define CHECKMAT(mat) \
	{bool flag = false; \
	for (int i = 0 ; i < mat.rows ; i ++) \
		for (int j = 0 ; j < mat.cols ; j ++) \
			if (std::isnan(mat( i, j))) flag = true;\
	if (flag) cout << "NaN : "<< #mat << endl; }

#define PRINTIF(mat1, mat2) 											\
	{ bool flag = false; 												\
	for (int i = 0 ; i < mat1.rows ; i ++)								\
		for (int j = 0 ; j < mat1.cols ; j ++) 						\
			if (std::isnan(mat1( i, j))) flag = true; 				\
	if (flag) {cout << "Mat Nan : " << #mat1 << " printMat : " << #mat2 << endl;	\
	for (int i = 0 ; i < mat2.rows ; i ++) 							\
		for (int j = 0 ; j < mat2.cols ; j ++) 						\
			cout << mat2( i, j) << "\t";							\
	cout << endl;}}



void FastHessian::FittingQuadric(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b,
                                  double* xX, double* xY, double* xZ, double* xS, bool corner )
{

	// Fitting
//	double x[4] = { 0 };
	cv::Matx< double, 4, 1 > X;

	cv::Matx< double, 81, 1 > Resp;
	cv::Matx< double, 81, 15 > M;

	int index  = 0;
	int index2 = 0;

	for (int i = -1 ; i <= 1 ; i++)
	for (int j = -1 ; j <= 1 ; j++)
	for (int k = -1 ; k <= 1 ; k++)
	{
		if (corner) {
			Resp( index2, 0 ) = b->getCornerResponse(r + i, c + j, d + k, t);
			index2++;
			Resp( index2, 0 ) = m->getCornerResponse(r + i, c + j, d + k, t);
			index2++;
			Resp( index2, 0 ) = t->getCornerResponse(r + i, c + j, d + k);
			index2++;
		} else {

		}

		for (int l = -1 ; l <= 1 ; l++)
		{
			// Matrice des Parametres
			M( index, 0 ) = i * i;
			M( index, 1 ) = j * j;
			M( index, 2 ) = k * k;
			M( index, 3 ) = l * l;
			M( index, 4 ) = i * j;
			M( index, 5 ) = i * k;
			M( index, 6 ) = i * l;
			M( index, 7 ) = j * k;
			M( index, 8 ) = j * l;
			M( index, 9 ) = k * l;
			M( index, 10 ) = i;
			M( index, 11 ) = j;
			M( index, 12 ) = k;
			M( index, 13 ) = l;
			M( index, 14 ) =  1;

			index++;
		}
	}

#define USE_SVD

#ifndef USE_SVD

	CvMat *M_t, *temp, *temp2, *temp3;

	cv::Matx<double 15, 81> M_t;
//	cvTranspose(M, M_t);
	M_t = M.t();

	cv::Matx< double, 15, 15 > temp;
//	cvMatMul(M_t, M, temp);
	temp = M_t * M

	cv::Matx< double, 15, 15 > temp2;
	temp2 = temp.inv();
//	cvInvert(temp, temp2);

	cv::Matx< doube, 15, 81> temp3;
//	cvMatMul(temp2, M_t, temp3);
	temp3 = temp2 * M_t;

	//cvInitMatHeader( &W, 15, 1, CV_64FC1, w, CV_AUTOSTEP );
	cv::Matx< double, 15, 1 > W;
	W = temp3 * Resp;
//	cvMatMul(temp3, Resp, &W);

#else
//{
//	CvMat *S, *S_inv, *S_inv_t, *V, *U_t, *tmp, *tmp2;

  cv::Matx< double, 81, 15> S;
  cv::Matx< double, 81, 15> S_inv;
  S_inv = cv::Matx< double, 81, 15>::zeros();
  cv::Matx< double, 15, 81> S_inv_t;
  cv::Matx< double, 15, 15> V;
  cv::Matx< double, 81, 81> U_t;

 cv::SVD::compute(M, S, U_t, V);

//  cvSVD(M, S, U_t, V, CV_SVD_U_T);

	CHECKMAT(M)
	CHECKMAT(S)
	CHECKMAT(U_t)
	CHECKMAT(V)
	CHECKMAT(Resp)


	bool flag = false;

	S_inv( 0, 0) =  1 / S( 0, 0 );

	for (int i = 1 ; i < 15 ; i++)
		if (S( i, i)/S( 0, 0) < 0.001)
			S_inv( i, i ) = 0;
		else
		{	S_inv( i, i ) = 1.0f/S( i, i );
			flag = true;}

	if (flag) cout << "truncated SVD" << endl;

	S_inv_t = S_inv.t();
//	cvTranspose(S_inv, S_inv_t);

	cv::Matx< double, 81, 1 > tmp;
//	tmp = cvCreateMat(81, 1, CV_64FC1);
	tmp = U_t * Resp;
//	cvMatMul(U_t, Resp, tmp);

	PRINTIF(tmp, U_t)
	PRINTIF(tmp, Resp)

	cv::Matx< double, 15, 1 > tmp2;
//	cvMatMul(S_inv_t, tmp, tmp2);
	tmp2 = S_inv_t * tmp;
	CHECKMAT(tmp2)
//	cvInitMatHeader( &W, 15, 1, CV_64FC1, w, CV_AUTOSTEP );
	cv::Matx< double, 15, 1 > W;
//	cvMatMul(V, tmp2, &W);
	W = V * tmp2;
//}
#endif


	// Maximum
//	A = cvCreateMat( 4, 4, CV_64FC1 );
//	b_mat = cvCreateMat( 4, 1, CV_64FC1 );

	cv::Matx< double, 4, 4> A;
	cv::Matx< double, 4, 1> b_mat;

	CHECKMAT(W)

	//Matrice A
	A( 0, 0 ) = W( 0 ) * 2;
	A( 0, 1 ) = W( 4 );
	A( 0, 2 ) = W( 5 );
	A( 0, 3 ) = W( 6 );
	A( 1, 0 ) = W( 4 );
	A( 1, 1 ) = W( 1 ) * 2;
	A( 1, 2 ) = W( 7 );
	A( 1, 3 ) = W( 9 );
	A( 2, 0 ) = W( 5 );
	A( 2, 1 ) = W( 8 );
	A( 2, 2 ) = W( 2 ) * 2;
	A( 2, 3 ) = W( 10 );
	A( 3, 0 ) = W( 6 );
	A( 3, 1 ) = W( 8 );
	A( 3, 2 ) = W( 9 );
	A( 3, 3 ) = W( 3 ) * 2;

	b_mat( 0 ) = W( 10 );
	b_mat( 1 ) = W( 11 );
	b_mat( 2 ) = W( 12 );
	b_mat( 3 ) = W( 13 );


#ifdef USE_SVD
{
//	CvMat *S, *S_inv, *V, *U_t, *tmp, *tmp2;

  cv::Matx< double, 4, 4 > S;
  cv::Matx< double, 4, 4 > S_inv = cv::Matx< double, 4, 4 >::zeros();
  cv::Matx< double, 4, 4 > V;
  cv::Matx< double, 4,4 > U_t;

  cv::SVD::compute( A, S, U_t, V);
//  cvSVD(A, S, U_t, V, CV_SVD_U_T);

	CHECKMAT(M)
	CHECKMAT(S)
	CHECKMAT(U_t)
	CHECKMAT(V)

	bool flag = false;

	S_inv( 0, 0 ) = 1.0 / S( 0, 0 );

	for (int i = 1 ; i < 4 ; i++)
		if ( S( i, i) /  S( 0, 0 ) < 0.001 )
			S_inv( i, i ) = 0;
		else
		{	S_inv( i, i ) = 1.0f / S( i, i );
			flag = true;}

	if (flag) cout << "truncated SVD 2" << endl;

	cv::Matx< double, 4, 1 > tmp;
//	cvMatMul(U_t, b_mat, tmp);
	tmp = U_t * b_mat;

	cv::Matx< double, 4, 1 > tmp2;
//	tmp2 = cvCreateMat(4, 1, CV_64FC1);
	tmp2 = S_inv * tmp;
//	cvMatMul(S_inv, tmp, tmp2);

//	cvInitMatHeader( &X, 4, 1, CV_64FC1, x, CV_AUTOSTEP );
	X = V * tmp2;
//	cvMatMul(V, tmp2, &X);

}

#else
	A_inv = cvCreateMat( 4, 4, CV_64FC1 );
	cvInvert(A, A_inv);
	cvInitMatHeader( &X, 4, 1, CV_64FC1, x, CV_AUTOSTEP );
	cvMatMul(A_inv, b_mat, &X);
#endif

	*xX = X( 0 );
	*xY = X( 1 );
	*xZ = X( 2 );
	*xS = X( 3 );

}


//-------------------------------------------------------

//! Computes the 3D Hessian matrix for a pixel.
cv::Matx44f* FastHessian::hessian4D(int d, int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, bool Corner)
{
	cv::Matx44f* H;
	double v, dxx, dyy, dzz, dss, dxy, dxz, dxs, dyz, dys, dzs;

	if (Corner) {
		v = m->getCornerResponse(r, c, d, t);

		dxx = m->getCornerResponse(r, c + 1, d, t) + m->getCornerResponse(r, c - 1, d, t) - 2 * v;
		dyy = m->getCornerResponse(r + 1, c, d, t) + m->getCornerResponse(r - 1, c, d, t) - 2 * v;
		dzz = m->getCornerResponse(r, c, d + 1, t) + m->getCornerResponse(r, c, d - 1, t) - 2 * v;
		dss = t->getCornerResponse(r, c, d) + b->getCornerResponse(r, c, d, t) - 2 * v;

		dxy = ( m->getCornerResponse(r + 1, c + 1, d, t) - m->getCornerResponse(r + 1, c - 1, d, t) -
				m->getCornerResponse(r - 1, c + 1, d, t) + m->getCornerResponse(r - 1, c - 1, d, t) ) / 4.0;
		dxz = ( m->getCornerResponse(r + 1, c, d + 1, t) - m->getCornerResponse(r + 1, c, d - 1, t) -
				m->getCornerResponse(r - 1, c, d + 1, t) + m->getCornerResponse(r - 1, c, d - 1, t) ) / 4.0;
		dxs = ( t->getCornerResponse(r, c + 1, d ) - t->getCornerResponse(r, c - 1, d) -
				b->getCornerResponse(r, c + 1, d, t) + b->getCornerResponse(r, c - 1, d, t) ) / 4.0;

		dyz = ( m->getCornerResponse(r, c + 1, d + 1, t) - m->getCornerResponse(r, c + 1, d - 1, t) -
				m->getCornerResponse(r, c - 1, d + 1, t) + m->getCornerResponse(r, c - 1, d - 1, t) ) / 4.0;
		dys = ( t->getCornerResponse(r + 1, c, d) - t->getCornerResponse(r - 1, c, d) -
				b->getCornerResponse(r + 1, c, d, t) + b->getCornerResponse(r - 1, c, d, t) ) / 4.0;

		dzs = ( t->getCornerResponse(r, c, d + 1) - t->getCornerResponse(r, c, d - 1) -
				b->getCornerResponse(r, c, d + 1, t) + b->getCornerResponse(r, c, d - 1, t) ) / 4.0;
	} else {
		v = m->getResponse(r, c, d, t);

		dxx = m->getResponse(r, c + 1, d, t) + m->getResponse(r, c - 1, d, t) - 2 * v;
		dyy = m->getResponse(r + 1, c, d, t) + m->getResponse(r - 1, c, d, t) - 2 * v;
		dzz = m->getResponse(r, c, d + 1, t) + m->getResponse(r, c, d - 1, t) - 2 * v;
		dss = t->getResponse(r, c, d) + b->getResponse(r, c, d, t) - 2 * v;

		dxy = ( m->getResponse(r + 1, c + 1, d, t) - m->getResponse(r + 1, c - 1, d, t) -
				m->getResponse(r - 1, c + 1, d, t) + m->getResponse(r - 1, c - 1, d, t) ) / 4.0;
		dxz = ( m->getResponse(r + 1, c, d + 1, t) - m->getResponse(r + 1, c, d - 1, t) -
				m->getResponse(r - 1, c, d + 1, t) + m->getResponse(r - 1, c, d - 1, t) ) / 4.0;
		dxs = ( t->getResponse(r, c + 1, d ) - t->getResponse(r, c - 1, d) -
				b->getResponse(r, c + 1, d, t) + b->getResponse(r, c - 1, d, t) ) / 4.0;

		dyz = ( m->getResponse(r, c + 1, d + 1, t) - m->getResponse(r, c + 1, d - 1, t) -
				m->getResponse(r, c - 1, d + 1, t) + m->getResponse(r, c - 1, d - 1, t) ) / 4.0;
		dys = ( t->getResponse(r + 1, c, d) - t->getResponse(r - 1, c, d) -
				b->getResponse(r + 1, c, d, t) + b->getResponse(r - 1, c, d, t) ) / 4.0;

		dzs = ( t->getResponse(r, c, d + 1) - t->getResponse(r, c, d - 1) -
				b->getResponse(r, c, d + 1, t) + b->getResponse(r, c, d - 1, t) ) / 4.0;
	}

  H = new cv::Matx44f;

  (*H)(0,0) = dxx;
  (*H)(0,1) = dxy;
  (*H)(0,2) = dxz;
  (*H)(0,3) = dxs;

  (*H)(1,0) = dxy;
  (*H)(1,1) = dyy;
  (*H)(1,2) = dyz;
  (*H)(1,3) = dys;

  (*H)(2,0) = dxz;
  (*H)(2,1) = dyz;
  (*H)(2,2) = dzz;
  (*H)(2,3) = dzs;

  (*H)(3,0) = dxs;
  (*H)(3,1) = dys;
  (*H)(3,2) = dzs;
  (*H)(3,3) = dss;

  return H;
}
/*
void FastHessian::print(CvMat* input)
{
	for (int x = 0 ; x < input->cols ; x++)
	{
		for (int y = 0 ; y < input->rows ; y++)
		{
			cout <<  (cvGet2D(input,y,x)).val[0] << "\t";

		}
		cout << endl;
	}

	cout << endl;
}
*/
void FastHessian::WriteResponseMap()
{

  for (std::vector<ResponseLayer *>::iterator it = responseMap.begin() ; it != responseMap.end(); ++it)
		{
			ResponseLayer * curr = *it;
			char filename[50];
			sprintf (filename, "FH/FastHessian_%d.", curr->filter);

			vtkImageData *Image = vtkImageData::New();
			Image->SetDimensions(curr->width, curr->height, curr->depth);
			Image->AllocateScalars(VTK_UNSIGNED_SHORT, 1);

			for (int i=0 ; i < curr->width ; i++)
				for (int j=0 ; j < curr->height ; j++)
					for (int k=0 ; k < curr->depth ; k++)
						*static_cast<unsigned short*>(Image->GetScalarPointer(i, j, k)) = 0;


			float* img_pointer = static_cast<float*>(Image->GetScalarPointer());


		  float max = VTK_FLOAT_MIN;
		  float min = VTK_FLOAT_MAX;

			int limit = 0; //std::max((float)ceil((float)((curr->filter - 1) / 2)/(float)curr->step)+2, (float)(0.6666*curr->filter + 2.5 * curr->step+2));

			for (int i=limit ; i < curr->width-limit ; i++)
				for (int j=limit ; j < curr->height-limit ; j++)
					for (int k=limit ; k < curr->depth-limit ; k++)
					{
						if (max < (curr->responses[i + j * curr->width + k * curr->width * curr->height]))
								max = (curr->responses[i + j * curr->width + k * curr->width * curr->height]);
						if (min > (curr->responses[i + j * curr->width + k * curr->width * curr->height]))
								min = (curr->responses[i + j * curr->width + k * curr->width * curr->height]);
					}

			float scale = (VTK_UNSIGNED_SHORT_MAX-VTK_UNSIGNED_SHORT_MIN-1)/(max-min);


			for (int i=limit ; i < curr->width-limit ; i++)
				for (int j=limit ; j < curr->height-limit ; j++)
					for (int k=limit ; k < curr->depth-limit ; k++)
						*static_cast<unsigned short*>(Image->GetScalarPointer(i, j, k)) = (unsigned short)
							((curr->responses[i + j * curr->width + k * curr->width * curr->height]- min)*scale);


			vtkMetaImageWriter *Writer=vtkMetaImageWriter::New();

			Image->Modified();

			Writer->SetInputData(Image);

			Writer->SetFileName(filename);

			Writer->Write();

			Writer->Delete();
			Image->Delete();
		}

}




//-------------------------------------------------------
