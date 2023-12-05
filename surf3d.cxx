#include <sstream>
#include <vector>
#include <map>
#include <vtkObjectFactory.h>
#include <vtkTimerLog.h>
#include <vtkMultiThreader.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageShiftScale.h>
#include <vtkImageResample.h>
#include <vtkImageLuminance.h>
#include <vtkImageData.h>
#include <vtkImageThreshold.h>
#include <vtkMetaImageReader.h>
#include <vtkNIFTIImageReader.h>
#include <vtkImageReader2Factory.h>
#include <string>
#include "integral.h"
#include "vtk3DSURF.h"
#include "vtkRobustImageReader.h"
#include "picojson.h"

using namespace std;


int main( int argc, char *argv[] )
{
	
	if (argc < 2) {
		cout << "Usage : surf3d file [options]" << endl;
		exit(1);
	}

	double spacing = 0;
	int maxSize = 0;
	double threshold = 0;
	char maskfilename[100];
	bool bmask = false;
	int writeJSON = 0;
	int writeBIN = 0;
	int writeCSV = 0;
	int writeCSVGZ = 1;
	int descriptorType = 0;
	int numberOfPoints = -1;
	int numberOfThreads = -1;
	int subVolumeRadius = 5;
	char *pointFile = 0;
	bool clampMinValues = false;
	float clampMinValue = 0;
	bool clampMaxValues = false;
	float clampMaxValue = 0;
	bool normalize = true;

	string outfilename("points");

	// Parse optionnal arguments
	int argumentsIndex = 2;
	while (argumentsIndex < argc) {
		char *key = argv[argumentsIndex];
		char *value = argv[argumentsIndex+1];

		if (strcmp(key,"-d") == 0) {
			maxSize = atoi(value);
		}
		if (strcmp(key,"-s") == 0) {
			spacing = atof(value);
		}

		if (strcmp(key,"-t") == 0) {
			threshold = atof(value);
		}

		if (strcmp(key,"-cmin") == 0) {
			clampMinValues = true;
			clampMinValue = atof(value);
		}

		if (strcmp(key,"-cmax") == 0) {
			clampMaxValues = true;
			clampMaxValue = atof(value);
		}

		if (strcmp(key,"-m") == 0) {
			strcpy(maskfilename, value);
			bmask = true;
		}

		if (strcmp(key,"-o") == 0) {
			outfilename = string(value );
		}

		if (strcmp(key,"-json") == 0) {
			writeJSON = atoi(value);
		}

		if (strcmp(key,"-csv") == 0) {
			writeCSV = atoi(value);
		}

		if (strcmp(key,"-bin") == 0) {
			writeBIN = atoi(value);
		}

		if (strcmp(key,"-csvgz") == 0) {
			writeCSVGZ = atoi(value);
		}

		if (strcmp(key,"-type") == 0) {
			descriptorType = atoi(value);
		}

		if (strcmp(key,"-n") == 0) {
			numberOfPoints = atoi(value);
		}

		if (strcmp(key,"-nt") == 0) {
			numberOfThreads = atoi(value);
		}

		if (strcmp(key,"-p") == 0) {
			pointFile = value;
		}

		if (strcmp(key,"-r") == 0) {
			subVolumeRadius = atoi(value);
		}

		if (strcmp(key,"-normalize") == 0) {
			normalize = atoi(value);
		}

		argumentsIndex += 2;
	}

	cout << "load : " << argv[1] << endl;
	
	vtkTimerLog *Timer = vtkTimerLog::New();
	Timer->StartTimer();

    vtkRobustImageReader *imageReader = vtkRobustImageReader::New();
	imageReader->SetFileName(argv[1]) ;
	imageReader->Update() ;
	vtkImageData *image = imageReader->GetOutput();
	vtkImageData *mask = 0;
	
	if (bmask)
	{
			// Create image reader factory and register Meta image reader with it
		vtkImageReader2Factory *maskReaderFactory = vtkImageReader2Factory::New();
		vtkMetaImageReader *metamaskReader = vtkMetaImageReader::New();
		maskReaderFactory->RegisterReader(metamaskReader);
		
		// Create a reader for image and try to load it
		vtkImageReader2* maskReader = maskReaderFactory->CreateImageReader2(maskfilename) ;
		if (!maskReader) {
			cerr << "Cannot load file " << argv[1] << " as an mask file; terminating.\n" ;
			return 5 ;
		}
		
		maskReader->SetFileName(maskfilename) ;
		maskReader->Update() ;
		mask = maskReader->GetOutput();
		
		int* mask_dims  = mask->GetDimensions();
		int* image_dims = image->GetDimensions();
		
		if (mask_dims[0] != image_dims[0] ||
				mask_dims[1] != image_dims[1] ||
				mask_dims[2] != image_dims[2] )
			{
				cerr << "Mismatch between Image && mask dimension" << endl;
				return 6;
			}
			
		if (mask->GetNumberOfScalarComponents() != 1)
			{
				cerr << "Too many components for mask (>1)" << endl;
				return 6;
			}
	}

	if ( clampMinValues ) {

		vtkImageThreshold *threshold = vtkImageThreshold::New();
		threshold->ReplaceOutOn();
		threshold->ThresholdByUpper( clampMinValue );
		threshold->SetOutValue( clampMinValue );
		threshold->SetInputData( image );
		threshold->Update();
		image = threshold->GetOutput();

	}

	if ( clampMaxValues ) {

		vtkImageThreshold *threshold = vtkImageThreshold::New();
		threshold->ReplaceOutOn();
		threshold->ThresholdByLower( clampMaxValue );
		threshold->SetOutValue( clampMaxValue );
		threshold->SetInputData( image );
		threshold->Update();
		image = threshold->GetOutput();
	}

	Timer->StopTimer();
	cout << "Image loaded in " << Timer->GetElapsedTime() << "s" << endl;
	
	vtk3DSURF *SURF = vtk3DSURF::New();
	SURF->SetInput(image);
	SURF->SetNbThread(24);
	SURF->SetNormalize(normalize);
	SURF->SetMaxSize(maxSize);
	SURF->SetSpacing(spacing);
	SURF->SetThreshold(threshold);
	SURF->SetDescriptorType(descriptorType);
	SURF->SetNumberOfPoints(numberOfPoints);
	SURF->SetSubVolumeRadius(subVolumeRadius);
	SURF->SetNbThread(numberOfThreads);

	double boundsArray[ 6 ];
	image->GetBounds( boundsArray );
	picojson::object root, bounds;
	bounds[ "xmin" ] = picojson::value( boundsArray[ 0 ] );
	bounds[ "xmax" ] = picojson::value( boundsArray[ 1 ] );
	bounds[ "ymin" ] = picojson::value( boundsArray[ 2 ] );
	bounds[ "ymax" ] = picojson::value( boundsArray[ 3 ] );
	bounds[ "zmin" ] = picojson::value( boundsArray[ 4 ] );
	bounds[ "zmax" ] = picojson::value( boundsArray[ 5 ] );
	root[ "bounds" ] = picojson::value( bounds );

	ofstream boundsFile;
	boundsFile.open( outfilename+".json" , std::ofstream::out | std::ofstream::trunc );
	boundsFile << picojson::value( root );
	boundsFile.close();

	if ( pointFile ) {

		cout << "Use points in " << pointFile << endl;
		SURF->SetPointFile( pointFile );

	}
	
	if (bmask) {

		SURF->SetMask(mask);
		cout << "Yeah, We have a mask !" << endl;

	}

	SURF->Update();

	if ( writeJSON ) {

		Timer->StartTimer();
		SURF->WritePoints( (outfilename+".json").c_str() );
		Timer->StopTimer();
		cout << "json written in " << Timer->GetElapsedTime() << "s" << endl;

	}

	if ( writeBIN ) {

		Timer->StartTimer();
		SURF->WritePointsBinary( (outfilename+".bin").c_str() );
		Timer->StopTimer();
		cout << "bin written in " << Timer->GetElapsedTime() << "s" << endl;

	}

	if ( writeCSV ) {

		Timer->StartTimer();
		SURF->WritePointsCSV( (outfilename+".csv").c_str() );
		Timer->StopTimer();
		cout << "csv written in " << Timer->GetElapsedTime() << "s" << endl;

	}

	if ( writeCSVGZ ) {

		Timer->StartTimer();
		SURF->WritePointsCSVGZ( (outfilename+".csv.gz").c_str() );
		Timer->StopTimer();
		cout << "csvgz written in " << Timer->GetElapsedTime() << "s" << endl;

	}

	return 0;
}
