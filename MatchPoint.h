#include <sstream>
#include <vector>
#include <map>
#include "ipoint.h"
#include <iostream>
#include <fstream>
#include "picojson.h"
#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <TooN/sim3.h>
#include <TooN/SVD.h>

using namespace picojson;
using namespace TooN;

typedef struct{
    IpVec points;
    int dimensions[3];
    double origin[3];
    double spacing[3];
} JSONdata;


class Box3 {
public:
	Box3();
	void AddPoint(TooN::Vector<3> point);
	bool ContainsPoint(Ipoint& p);
	TooN::Vector<3> min;
	TooN::Vector<3> max;
};

class MatchPoint {

public:
  //! Constructor
  MatchPoint():RansacDist(50), 
				RansacMinInliers(2), 
				MatchingDist(0.95), 
				MatchingDist2Second(0.98),
				MatchingScale(1.5),
				echec(false),
				useBBoxin(false) {};
  
  
	    //! Save the parameters
    void saveParameters(const float aRansacDist, 
                        const int 	aRansacMinInliers,
                        const float aMatchingDist, 
                        const float aMatchingDist2Second,
                        const float aMatchingScale);
                        
	void Parse(const char *fileName, int id);
	
	void computeMatches();
	
	void computeTransform();
	
	void BboxParse(std::string bbox);
	
	
	void WriteTransform(const char *fileName, bool writeInliers );

	void setComputeBoundingBoxes(bool value) {
		this->computeBoundingBoxes = value;
	}

private:
	bool getRT(std::map<int, std::pair<Ipoint, Ipoint>* > Pairs, TooN::SIM3<double>& Transform);

	JSONdata volumes[2];
	float RansacDist;
	int RansacMinInliers;
	float MatchingDist;
	float MatchingDist2Second;
	float MatchingScale;
	IpPairVec matches;
	SIM3<double> BestTransform;
	int maxinlier;
	bool echec;
	Box3 bboxAmax, bboxBmax;
	Box3 bboxAin, bboxBin;
	int nbPointInA, nbPointInB;
	std::vector< int > inliers;
	bool useBBoxin;
	bool computeBoundingBoxes;
};

std::string base64_decode(std::string const& s);

