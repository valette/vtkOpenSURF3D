#include <sstream>
#include "MatchPoint.h"

using std::cout, std::endl;

int main( int argc, char *argv[] )
{
	if (argc < 3) {
		cout << "Usage : SURF file1 file2 [options]" << endl;
		exit(1);
	}
	
	int RansacMinInliers(2); 

	double RansacDist(50), 
			MatchingDist(0.95), 
			MatchingDist2Second(0.98),
			MatchingScale(1.5);

	bool computeBoxes = false;
	bool writeInliers = false;

	// Parse optionnal arguments
	int argumentsIndex = 3;
	std::string bbox;
	
	while (argumentsIndex < argc) {
		char* key = argv[argumentsIndex];
		char *value = argv[argumentsIndex + 1];

		if (strcmp(key, "-Rd") == 0) {
			RansacDist = atof(value);
		}
		if (strcmp(key, "-Ri") == 0) {
			RansacMinInliers = atoi(value);
		}
		if (strcmp(key, "-Md") == 0) {
			MatchingDist = atof(value);
		}
		if (strcmp(key, "-Md2") == 0) {
			MatchingDist2Second = atof(value);
		}
		if (strcmp(key, "-Ms") == 0) {
			MatchingScale = atof(value);
		}
		if (strcmp(key, "-b") == 0) {
			computeBoxes = atoi(value);
		}
		if (strcmp(key, "-bb") == 0) {
			bbox = value;
		}
		if (strcmp(key, "-i") == 0) {
			writeInliers = true;
			argumentsIndex -= 1;
		}
		
		argumentsIndex += 2;
	}
	MatchPoint mp;
	mp.setComputeBoundingBoxes(computeBoxes);
	cout << "Save Parameters" << endl;
	mp.saveParameters(RansacDist, RansacMinInliers, MatchingDist, MatchingDist2Second, MatchingScale);
	cout << "Parse " << argv[1] << endl;
	mp.Parse(argv[1], 0);
	cout << "Parse " << argv[2] << endl;
	mp.Parse(argv[2], 1);	
	
	if (!bbox.empty()) {
		cout << "Parse Bounding Box" << endl;
		mp.BboxParse(base64_decode(bbox));
	}
	
	cout << "compute Matches" << endl;
	mp.computeMatches();
	cout << "compute Transform" << endl;
	mp.computeTransform();
	cout << "write JSON file" << endl;
	mp.WriteTransform("transform.json", writeInliers );
	cout<<endl;
}
