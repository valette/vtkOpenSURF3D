#include <string>
#include <vtkTimerLog.h>
#include <vtkImageData.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageMirrorPad.h>
#include <vtkImageThreshold.h>
#include <vtkImageTranslateExtent.h>
#include <vtkMetaImageReader.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>

#include "vtk3DSURF.h"
#include "vtkRobustImageReader.h"
#include "picojson.h"

int main( int argc, char *argv[] )
{
	
	if (argc < 2) {
		cout << "Usage : surf3d file [options]" << endl;
		cout << "Available options:" << endl;
		cout << "-bin 0/1       : write points as bin file. Default : 0" << endl;
		cout << "-csv 0/1       : write points as csv file. Default : 0" << endl;
		cout << "-cmin value    : clamp values lower than specified value" << endl;
		cout << "-cmax value    : clamp values larger than specified value" << endl;
		cout << "-csvgz 0/1     : write points as csv.gz file. Default : 1" << endl;
		cout << "-gz opts       : set gz options such as compression level" << endl;
		cout << "-precision n   : set coefficients precision in csv.gz file" << endl;
		cout << "-json 0/1      : write points as json file. Default : 0" << endl;
		cout << "-n number      : maximum number of points" << endl;
		cout << "-normalize 0/1 : normalize descriptors (default : 1 )" << endl;
		cout << "-o basename    : set output file name. Default: \"points\"" << endl;
		cout << "-r radius      : descriptor volume radius. Default : 5" << endl;
		cout << "-s spacing     : resample input image to isotropic sampling with given spacing" << endl;
		cout << "-t threshold   : set detector threshold. Default: 0" << endl;
		cout << "-type 0/1/2    : set descriptor type :" << endl;
		cout << "       0 : SURF3D descriptor (default). Descriptor size : 48" << endl;
		cout << "       1 : subvolume HAAR coefficients. Descriptor size : 24 * radius^3" << endl;
		cout << "       2 : subvolume raw voxels. Descriptor size : 8 * radius^3" << endl;
		exit(1);
	}

	double spacing = 0;
	int maxSize = 0;
	double threshold = 0;
	char *maskfilename = 0;
	int writeJSON = 0;
	int writeBIN = 0;
	int writeCSV = 0;
	int writeCSVGZ = 1;
	int descriptorType = 0;
	int numberOfPoints = -1;
	int numberOfThreads = -1;
	int pad = 0;
	int subVolumeRadius = 5;
	char *pointFile = 0;
	bool clampMinValues = false;
	float clampMinValue = 0;
	bool clampMaxValues = false;
	float clampMaxValue = 0;
	bool normalize = true;
	char *gzOpts = 0;
	int precision = -1;

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
			maskfilename = value;
		}

		if (strcmp(key,"-o") == 0) {
			outfilename = value;
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

		if (strcmp(key,"-pad") == 0) {
			pad = atoi(value);
		}

		if (strcmp(key,"-gz") == 0) {
			gzOpts = value;
		}

		if (strcmp(key,"-precision") == 0) {
			precision = atoi( value );
		}

		argumentsIndex += 2;
	}

	cout << "load : " << argv[1] << endl;
	
	vtkNew<vtkTimerLog> Timer;
	Timer->StartTimer();

    vtkNew<vtkRobustImageReader> imageReader;
	imageReader->SetFileName(argv[1]) ;
	imageReader->Update() ;
	vtkSmartPointer<vtkImageData> image = imageReader->GetOutput();
	vtkSmartPointer<vtkImageData> mask;

	if (maskfilename)
	{
			// Create image reader factory and register Meta image reader with it
		vtkNew<vtkImageReader2Factory> maskReaderFactory;
		vtkNew<vtkMetaImageReader> metamaskReader;
		maskReaderFactory->RegisterReader(metamaskReader);
		
		// Create a reader for image and try to load it
		vtkSmartPointer<vtkImageReader2> maskReader = maskReaderFactory->CreateImageReader2(maskfilename) ;
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

	if ( pad ) {
		std::cout << "Apply mirror padding, with padding = " << pad << std::endl;
		int dimensions[ 3 ];
		double spacing[ 3 ];
		image->GetDimensions( dimensions );
		image->GetSpacing( spacing );
		vtkNew<vtkImageMirrorPad> padding;
		padding->SetInputData( image );
		int p[ 3 ];
		for ( int i = 0; i < 3; i++ ) p[ i ] = ceil( pad / spacing[ i ] );
		std::cout << "Voxel padding = " << p[ 0 ] << " " << p[ 1 ] << " " << p[ 2 ] << std::endl;
		padding->SetOutputWholeExtent( -p[ 0 ], dimensions[ 0 ] + p[ 0 ],
			-p[ 1 ], dimensions[ 1 ] + p[ 1 ], -p[ 2 ], dimensions[ 2 ] + p[ 2 ] );
		padding->Update();
		vtkNew<vtkImageTranslateExtent> translate;
		translate->SetInputData( padding->GetOutput() );
		translate->SetTranslation( p[ 0 ], p[ 1 ], p[ 2 ] );
		translate->Update();
		image = translate->GetOutput();
	}

	if ( clampMinValues ) {

		vtkNew<vtkImageThreshold> threshold;
		threshold->ReplaceOutOn();
		threshold->ThresholdByUpper( clampMinValue );
		threshold->SetOutValue( clampMinValue );
		threshold->SetInputData( image );
		threshold->Update();
		image = threshold->GetOutput();

	}

	if ( clampMaxValues ) {

		vtkNew<vtkImageThreshold> threshold;
		threshold->ReplaceOutOn();
		threshold->ThresholdByLower( clampMaxValue );
		threshold->SetOutValue( clampMaxValue );
		threshold->SetInputData( image );
		threshold->Update();
		image = threshold->GetOutput();
	}

	Timer->StopTimer();
	cout << "Image loaded in " << Timer->GetElapsedTime() << "s" << endl;
	
	vtkNew<vtk3DSURF> SURF;
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
	
	if (maskfilename) {

		cout << "Use mask : " << maskfilename << endl;
		SURF->SetMask(mask);

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
		SURF->WritePointsCSVGZ( (outfilename+".csv.gz").c_str(), gzOpts, precision );
		Timer->StopTimer();
		cout << "csvgz written in " << Timer->GetElapsedTime() << "s" << endl;

	}

	return 0;
}
