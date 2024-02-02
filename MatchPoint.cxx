#include <sstream>
#include <vector>
#include <map>
#include "ipoint.h"
#include <iostream>
#include <fstream>
#include "picojson.h"
#include "MatchPoint.h"
#include <float.h>
#include <stdlib.h>     /* srand, rand */

using namespace picojson;
using namespace std;

void MatchPoint::Parse(const char *fileName, int id) {
	
	if (id != 0 && id != 1)
	{
	  std::cerr << "id must be 0 or 1" << std::endl;
	  exit(2);
	}

	ifstream pointsFile(fileName, ios::in);
	picojson::value v;
	pointsFile >> v;
	
	// parse the input
	std::string err = picojson::get_last_error();
	if (! err.empty()) {
	  std::cerr << err << std::endl;
	  exit(1);
	}

	// check if the type of the value is "object"
	if (! v.is<picojson::object>()) {
	  std::cerr << "JSON is not an object" << std::endl;
	  exit(2);
	}

	// obtain a const reference to the map, and print the contents
    picojson::object& o = v.get<picojson::object>();
    
	for (const auto &j : o["points"].get<picojson::array>()) {
		Ipoint point;
		picojson::object p = j.get<picojson::object>();
		point.laplacian = p["laplacian"].get<double>();
		point.response  = p["response"].get<double>();
		point.scale     = p["scale"].get<double>();
		point.x 		= p["x"].get<double>();
		point.y 		= p["y"].get<double>();
		point.z 		= p["z"].get<double>();
		picojson::array t = p["descriptor"].get<picojson::array>();
		point.allocate( t.size() );
		for ( const auto & m : t ) point.descriptor.push_back( m.get<double>() );
		this->points[id].push_back(point);
    }

}

//! Populate IpPairVec with matched ipts 
void MatchPoint::computeMatches()
{
	IpVec &ipts1 = points[0];
	IpVec &ipts2 =  points[1];
	
	cout << "Set1 size : " << ipts1.size() << endl;
	cout << "Set2 size : " << ipts2.size() << endl;

	if (useBBoxin)
	{
		int size = 0;
		for( const auto &pt : ipts1 )
			if ( bboxAin.ContainsPoint( pt ) ) size++;

		cout << "Set1 size (in bbox) : " << size << endl;
		size = 0;
		for( const auto &pt : ipts2 )
			if ( bboxBin.ContainsPoint( pt ) )size++;
		cout << "Set2 size (in bbox) : " << size << endl;
	}

    float dist, d1, d2;
    int bestJ;
    matches.clear();
    if ((ipts1.size() < 2) || (ipts2.size() < 2)) return;

	cout << "Config : " << endl;
	cout << "\tRansacDist :: " << RansacDist << endl;
	cout << "\tRansacMinInliers :: " << RansacMinInliers << endl;
	cout << "\tMatchingDist :: " << MatchingDist << endl;
	cout << "\tMatchingDist2Second :: " << MatchingDist2Second << endl;
	cout << "\tMatchingScale :: " << MatchingScale << endl;

  for(unsigned int i = 0; i < ipts1.size(); i++) 
  {
	  //cout << ipts1[i].x << "\t" << ipts1[i].y << "\t" << ipts1[i].z << endl;
	  
	  if (!bboxAin.ContainsPoint(ipts1[i]) && useBBoxin)
		continue;
	  
	  
    d1 = d2 = FLT_MAX;


    for(unsigned int j = 0; j < ipts2.size(); j++) 
    {
	  if (!bboxBin.ContainsPoint(ipts2[j]) && useBBoxin)
		continue;
		
      if (ipts1[i].laplacian != ipts2[j].laplacian)
		  continue;
	  
	  
	if ((ipts1[i].scale/ipts2[j].scale > MatchingScale) || 
			(ipts2[j].scale/ipts1[i].scale > MatchingScale) )
			continue;

	  
      dist = ipts1[i] - ipts2[j];

      if(dist<d1) // if this feature matches better than current best
      {
        d2 = d1;
        d1 = dist;
        bestJ = j;
      }
      else if(dist<d2) // this feature matches better than second best
      {
        d2 = dist;
      }
    }

    if (	(d1/d2 < MatchingDist2Second) && (d2 != FLT_MAX) && 
			(d1 < MatchingDist))
    { 
      matches.push_back(std::make_pair(i, bestJ ));
    }
  }
    cout << "Nb match : " << matches.size() << endl;
}

void MatchPoint::saveParameters(const float aRansacDist, 
                        const int 	aRansacMinInliers,
                        const float aMatchingDist, 
                        const float aMatchingDist2Second,
                        const float aMatchingScale) {
			RansacDist	= aRansacDist;
			RansacMinInliers = aRansacMinInliers;
			MatchingDist = aMatchingDist;
			MatchingDist2Second = aMatchingDist2Second;
			MatchingScale = aMatchingScale;
}



void MatchPoint::computeTransform()
{
	int N = 3;
	const auto &pts1 = this->points[ 0 ];
	const auto &pts2 = this->points[ 1 ];

	if (matches.size() <= N) {
		cout << "FAIL 1 : " << matches.size() << endl;
		echec = true;
		return;
	}

	int fail = 0;
	TooN::SIM3<double> Transform;
	maxinlier = 0;
	srand (1000);
	int *bestInliers = new int[ 2 * matches.size() ];
	int *currentInliers = new int[ 2 * matches.size() ];
	
	for (int i = 0; i < 8000 && fail < 10000 ; i++)
	{
		fail++;

		// - 1 - Pick up 2 points and compute transformation
		std::set<int> RANSACSet;

		while (RANSACSet.size() < N) {
			int x = rand() % matches.size();
			if (RANSACSet.count(x) > 0)	continue;
			RANSACSet.insert( x );
		}

		TooN::SIM3<double> Transform;

		if (!getRT(RANSACSet, Transform))
			{
				cout << "Fail" << endl;
				i--;
				continue;
			}
		/*else
			{
				cout << Transform << endl;
			}*/


		int inliers = 0;
		Box3 bboxA, bboxB;
		TooN::SIM3<double> TransformInv = Transform.inverse();
		
		TooN::Vector<3> pointA;
		TooN::Vector<3> pointB;

		// - 2 - Calcul inliers
		for (auto it : matches )
		{
			//pointA = makeVector(it->first.x, it->first.y, it->first.z);
			const Ipoint &ptA = pts1[ it.first ];
			const Ipoint &ptB = pts2[ it.second ];
			pointA[0] = ptA.x;
			pointA[1] = ptA.y;
			pointA[2] = ptA.z;
			
			//pointB = makeVector(it->second.x, it->second.y, it->second.z);
			pointB[0] = ptB.x;
			pointB[1] = ptB.y;
			pointB[2] = ptB.z;
			
			//cout << TooN::norm( pointB - Transform*pointA ) << endl;
			
			if (TooN::norm( pointB - TransformInv*pointA ) < RansacDist) {
				currentInliers[ 2 * inliers ] = it.first;
				currentInliers[ 1 + 2 * inliers ] = it.second;
				inliers++;
				if (this->computeBoundingBoxes) {
					bboxA.AddPoint(pointA);
					bboxB.AddPoint(pointB);
				}
			}
		}
		
		if (inliers > maxinlier)
		{
			maxinlier = inliers;
			BestTransform = TransformInv;
			int *temp = currentInliers;
			currentInliers = bestInliers;
			bestInliers = temp;
			bboxAmax = bboxA;
			bboxBmax = bboxB;
		}
	}

	nbPointInA = 0;
	nbPointInB = 0;

	for ( auto it : matches ) {

		if ( bboxAmax.ContainsPoint( pts1[ it.first ] ) ) nbPointInA++;
		if ( bboxBmax.ContainsPoint( pts2[ it.second ] ) ) nbPointInB++;

	}

	if ( maxinlier < RansacMinInliers ) echec = true;
	cout << "Max inliers : " << maxinlier << endl;
	this->inliers.clear();
	for ( int i = 0; i < 2 * maxinlier; i++ )
		this->inliers.push_back( bestInliers[ i ] );

	delete [] bestInliers;
	delete [] currentInliers;
}


bool MatchPoint::getRT(const std::set<int> &Pairs, TooN::SIM3<double>& Transform)
{
    //[Least-Squares Estimation... Umeyama]
	const int nbrFinalInliers = Pairs.size();

	TooN::Matrix<Dynamic, 3> A(nbrFinalInliers, 3);
	TooN::Matrix<Dynamic, 3> B(nbrFinalInliers, 3);

	//remplir A et B

	TooN::Vector<3> moyA = Zeros;
	TooN::Vector<3> moyB = Zeros;

	unsigned int row = 0;

	for ( const auto &index : Pairs ) {
		const auto &pair = this->matches[ index ];
		const auto &ptA = this->points[ 0 ][ pair.first ];
		const auto &ptB = this->points[ 1 ][ pair.second ];
		TooN::Vector<3> pointA = makeVector( ptA.x, ptA.y, ptA.z);
		TooN::Vector<3> pointB = makeVector( ptB.x, ptB.y, ptB.z);
		
		A[row] = pointA;
		B[row] = pointB;

		moyA += pointA;
		moyB += pointB;
		row++;
	}

	moyA /= nbrFinalInliers;
	moyB /= nbrFinalInliers;
	
	int nb = row;
	
	for (row = 0 ; row < nb ; row++)
	{
		A[row] -= moyA;
		B[row] -= moyB;
	}
	
	TooN::Matrix<3> Ca = A.T()*A;
	TooN::Matrix<3> Cb = B.T()*B;
	
	TooN::SVD<3> svdCa(Ca);
	TooN::SVD<3> svdCb(Cb);
	
	TooN::Vector<3> EigCa = svdCa.get_diagonal();
	TooN::Vector<3> EigCb = svdCb.get_diagonal();
	
	double s = sqrt( (EigCa*EigCb)/(EigCa*EigCa) );
	
/*
	TooN::Matrix<3> C = A.T()*B;

	TooN::SVD<3> svdC(C);

	TooN::Matrix<3> VT = svdC.get_VT();
	TooN::Matrix<3> U = svdC.get_U();

	double detuv = TooN::determinant(U)*TooN::determinant(VT); //cv::determinant(svdC.u)*cv::determinant(svdC.vt);

	TooN::Matrix<3> R = Identity;

	if (detuv >0)
	{
		R = U*VT;
	}
	else
	{
		return false;
	}
*/
	TooN::Matrix<3> R = Identity;
	
    Transform.get_translation() = moyA - R*moyB;
    Transform.get_rotation() 	= R;
    Transform.get_scale() 		= s;

	return true;
}

void MatchPoint::WriteTransform(const char *fileName, bool writeInliers) {

	picojson::object bboxA;
	picojson::object bboxB;
	
	picojson::array abboxAmin;
	picojson::array abboxAmax;
	picojson::array abboxBmin;
	picojson::array abboxBmax;
	
	picojson::array translation;
	picojson::array rotation;
	TooN::Vector<> vector = BestTransform.get_translation();
	TooN::Matrix<> matrice = BestTransform.get_rotation().get_matrix();
	double scale = BestTransform.get_scale();
	
	for (int i = 0; i < 3; i++) {
		abboxAmin.push_back(value(bboxAmax.min[i]));
		abboxAmax.push_back(value(bboxAmax.max[i]));
		abboxBmin.push_back(value(bboxBmax.min[i]));
		abboxBmax.push_back(value(bboxBmax.max[i]));
		
		translation.push_back(value(vector[i]));
	}
		
	bboxA["min"] = value(abboxAmin);
	bboxA["max"] = value(abboxAmax);
	bboxB["min"] = value(abboxBmin);
	bboxB["max"] = value(abboxBmax);
	
	for (int i = 0; i < 3; i++)
	{
		picojson::array row;
		for (int j = 0; j < 3; j++)
			row.push_back(value(matrice[i][j]));
		rotation.push_back(value(row));
	}
		
	object root;
	root["translation"] = value(translation);
	root["rotation"] 	= value(rotation);
	root["scale"]		= value(scale);
	root["inliers"]		= value((float)maxinlier);
	root["fail"]		= value(echec);
	root["bboxA"]		= value(bboxA);
	root["bboxB"]		= value(bboxB);
	root["nbPointInA"]  = value((float)nbPointInA);
	root["nbPointInB"]  = value((float)nbPointInB);

	if ( writeInliers ) {

		picojson::array allInliers;

		for ( int i = 0; i < this->inliers.size() / 2; i++ ) {

			picojson::array pair;
			pair.push_back( value( ( double ) this->inliers[ 2 * i ] ) );
			pair.push_back( value( ( double ) this->inliers[ 1 + 2 * i ] ) );
			allInliers.push_back( value( pair ) );

		}

		root["allInliers"]  = value( allInliers );

	}

	ofstream pointsFile;
	pointsFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	pointsFile << value(root);
	std::cerr << value(root);
	pointsFile.close();
}

void MatchPoint::BboxParse(string bbox) {
	picojson::value v;
	
	const char* p = bbox.data();
	
	cout << p << endl;
	picojson::parse(v, p, p + strlen(p));
	
	// check if the type of the value is "object"
	if (! v.is<picojson::object>()) {
	  std::cerr << "JSON bbox is not an object" << std::endl;
	  exit(2);
	}
	
    picojson::object& o = v.get<picojson::object>();
    
	bboxAin.min[0] = o["bboxA"]
						.get<picojson::object>()["min"]
						.get<picojson::object>()["x"]
						.get<double>();
	bboxAin.min[1] = o["bboxA"]
						.get<picojson::object>()["min"]
						.get<picojson::object>()["y"]
						.get<double>();
	bboxAin.min[2] = o["bboxA"]
						.get<picojson::object>()["min"]
						.get<picojson::object>()["z"]
						.get<double>();

	bboxBin.min[0] = o["bboxB"]
						.get<picojson::object>()["min"]
						.get<picojson::object>()["x"]
						.get<double>();
	bboxBin.min[1] = o["bboxB"]
						.get<picojson::object>()["min"]
						.get<picojson::object>()["y"]
						.get<double>();
	bboxBin.min[2] = o["bboxB"]
						.get<picojson::object>()["min"]
						.get<picojson::object>()["z"]
						.get<double>();
	
	bboxAin.max[0] = o["bboxA"]
						.get<picojson::object>()["max"]
						.get<picojson::object>()["x"]
						.get<double>();
	bboxAin.max[1] = o["bboxA"]
						.get<picojson::object>()["max"]
						.get<picojson::object>()["y"]
						.get<double>();
	bboxAin.max[2] = o["bboxA"]
						.get<picojson::object>()["max"]
						.get<picojson::object>()["z"]
						.get<double>();

	bboxBin.max[0] = o["bboxB"]
						.get<picojson::object>()["max"]
						.get<picojson::object>()["x"]
						.get<double>();
	bboxBin.max[1] = o["bboxB"]
						.get<picojson::object>()["max"]
						.get<picojson::object>()["y"]
						.get<double>();
	bboxBin.max[2] = o["bboxB"]
						.get<picojson::object>()["max"]
						.get<picojson::object>()["z"]
						.get<double>();
						
	useBBoxin = true;
	computeBoundingBoxes = false;
}


Box3::Box3() {
	min = makeVector(1e10, 1e10, 1e10);
	max = makeVector(-1e10, -1e10, -1e10);
	
}

void Box3::AddPoint(const TooN::Vector<3> &point) {
	min[0] = std::min(point[0], min[0]);
	min[1] = std::min(point[1], min[1]);
	min[2] = std::min(point[2], min[2]);
	
	max[0] = std::max(point[0], max[0]);
	max[1] = std::max(point[1], max[1]);
	max[2] = std::max(point[2], max[2]);
}

bool Box3::ContainsPoint(const Ipoint& p) {
		if ( p.x < min[0] || p.x > max[0] ||
		     p.y < min[1] || p.y > max[1] ||
		     p.z < min[2] || p.z > max[2] ) {
				return false;
		}
	return true;
	
}

static const std::string base64_chars = 
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
             "abcdefghijklmnopqrstuvwxyz"
             "0123456789+/";


static inline bool is_base64(unsigned char c) {
  return (isalnum(c) || (c == '+') || (c == '/'));
}

std::string base64_decode(std::string const& encoded_string) {
  int in_len = encoded_string.size();
  int i = 0;
  int j = 0;
  int in_ = 0;
  unsigned char char_array_4[4], char_array_3[3];
  std::string ret;

  while (in_len-- && ( encoded_string[in_] != '=') && is_base64(encoded_string[in_])) {
    char_array_4[i++] = encoded_string[in_]; in_++;
    if (i ==4) {
      for (i = 0; i <4; i++)
        char_array_4[i] = base64_chars.find(char_array_4[i]);

      char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
      char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
      char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

      for (i = 0; (i < 3); i++)
        ret += char_array_3[i];
      i = 0;
    }
  }

  if (i) {
    for (j = i; j <4; j++)
      char_array_4[j] = 0;

    for (j = 0; j <4; j++)
      char_array_4[j] = base64_chars.find(char_array_4[j]);

    char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
    char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
    char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

    for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
  }

  return ret;
}
