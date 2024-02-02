#include <sstream>
#include <vector>
#include <map>
#include <vtkBoundingBox.h>
#include <vtkObjectFactory.h>
#include <vtkTimerLog.h>
#include <vtkMultiThreader.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageCast.h>
#include <vtkImageShiftScale.h>
#include <vtkImageResample.h>
#include <vtkImageLuminance.h>
#include <vtkImageResize.h>
#include <vtkMetaImageWriter.h>
#include <vtkNew.h>
#include "vtk3DSURF.h"
#include "integral.h"
#include <vtkVersion.h>

#include "fasthessian.h"
#include "ipoint.h"
#include "surf.h"

#include <zlib.h>

#include "picojson.h"

using namespace picojson;

vtkStandardNewMacro(vtk3DSURF);

bool compareResponses (Ipoint &i,Ipoint &j) { return (i.response>j.response); }

void vtk3DSURF::ReadIPoints() {

    std::string line;
    this->points.clear();
	std::cout << "Read : " << this->PointFile << std::endl;
	ifstream file( this->PointFile );

	double origin[ 3 ], spacing[ 3 ], bounds[ 6 ];
	vtkBoundingBox box;
	this->Cast->GetOrigin( origin );
	this->Cast->GetSpacing( spacing );
	box.SetBounds( this->Cast->GetBounds() );
	double Sspacing = pow( spacing[ 0 ] * spacing[ 1 ] * spacing[ 2 ], 1.0 / 3.0 );

    while( std::getline( file, line ) ) {

		std::stringstream lineStream(line);
		std::string  cell;

		Ipoint point;
		if ( !std::getline( lineStream, cell, ',' ) ) continue;
		float x = std::stof( cell );
		point.x = ( x - origin[ 0 ] ) / spacing[ 0 ];
		std::getline( lineStream, cell, ',' );
		float y = std::stof( cell );
		point.y = ( y - origin[ 1 ] ) / spacing[ 1 ];
		std::getline( lineStream, cell, ',' );
		float z = std::stof( cell );
		point.z = ( z - origin[ 2 ] ) / spacing[ 2 ];
		std::getline( lineStream, cell, ',' );
		point.scale = std::stof( cell ) / Sspacing;
		if ( !box.ContainsPoint( x, y, z ) ) {
			std::cout << "Error : point " << this->points.size();
			std::cout << " with coordinates : " << point.x << " " << point.y << " " << point.z << " is outside image" << std::endl;
		}
		this->points.push_back( point );

    }

}

void vtk3DSURF::Update() {

	vtkNew<vtkTimerLog> Timer;

	int dimensions[3];
	Input->GetDimensions(dimensions);
	vtkSmartPointer<vtkImageData> temp = Input;

	cout << "Image dimensions : " << dimensions[0] << " " <<
			dimensions[1] << " " << dimensions[2] << endl;

	if (temp->GetNumberOfScalarComponents() >= 3) {
		cout << "Image has several components, converting..." << endl;
		vtkNew<vtkImageLuminance> luminance;
		luminance->SetInputData(temp);
		luminance->Update();
		temp = luminance->GetOutput();
	}

	double originspacing[3];
	temp->GetSpacing(originspacing);
	cout << "Initial spacing  : " << originspacing[0] << " " << originspacing[1] << " " << originspacing[2] << endl;

	double newSpacing;
	if (Spacing != 0 || ((MaxSize > 0) && (dimensions[2] != 1)) ) {

		if (Spacing != 0) {

			if (MaxSize != 0)
				cerr << "/!\\ Do Not set Size and Spacing at the same time /!\\ " << endl;

			cout << "Final   spacing  : " << Spacing << endl;
			newSpacing = Spacing;

		} else if ((MaxSize > 0) && (dimensions[2] != 1)) {

			double spacing[3];
			temp->GetSpacing(spacing);
			double lmin = (double) dimensions[0] * spacing[0];
			int imin = 0, imin2 = 0;
			double lmin2 = lmin;

			for (int i = 1; i < 3; i++) {
				double length = (double) dimensions[i] * spacing[i];
				if (length < lmin) {
					lmin2 = lmin;
					imin2 = imin;

					lmin = length;
					imin = i;
				}
			}

			newSpacing = (double) lmin2 / MaxSize;
		}

		vtkNew<vtkImageResample> resample;
		resample->SetInputData(temp);
		for (int i = 0; i != 3; i++) {
			resample->SetAxisOutputSpacing(i, newSpacing);
		}
		if (this->NbThread > 0) {
			resample->SetNumberOfThreads(NbThread);
		}
		resample->Update();
		temp = resample->GetOutput();
		cout << "Image resampled with spacing = " << newSpacing << endl;
		temp->GetDimensions(dimensions);
		cout << "New image dimensions : " << dimensions[0] << " " <<
			dimensions[1] << " " << dimensions[2] << endl;

		if (Mask != 0)
		{
			vtkImageData* temp2;
			vtkNew<vtkImageResample> resamplemask;
			resamplemask->SetInputData(Mask);
			resamplemask->SetInterpolationModeToNearestNeighbor();
			for (int i = 0; i != 3; i++)
				resamplemask->SetAxisOutputSpacing(i, newSpacing);
			if (this->NbThread > 0)
				resamplemask->SetNumberOfThreads(NbThread);
			resamplemask->Update();
			this->Mask = resamplemask->GetOutput();
		}

	}

	vtkNew<vtkImageCast> imageCast;
	imageCast->SetClampOverflow(1);
	imageCast->SetInputData(temp);
	imageCast->SetOutputScalarTypeToInt ();
	if (this->NbThread > 0) imageCast->SetNumberOfThreads(NbThread);
	imageCast->Update();
	this->Resized = imageCast->GetOutput();

	double Range[2];
	temp->GetScalarRange(Range);
	cout << "Range " << Range[0] << " " <<  Range[1] << " " << " " << endl;
	vtkNew<vtkImageShiftScale> Shift;
	Shift->SetClampOverflow(1);
	Shift->SetInputData( this->Resized );
	Shift->SetOutputScalarTypeToInt ();
	Shift->SetShift( -Range[ 0 ] );
	if (this->NbThread > 0) Shift->SetNumberOfThreads(NbThread);
	Shift->Update();
	this->Cast = Shift->GetOutput();

	Timer->StartTimer();
	this->Integral =  vtkSmartPointer<vtkImageData>::Take(ComputeIntegral(Cast));
	Timer->StopTimer();
	cout << "Integral computed in " << Timer->GetElapsedTime() << "s" << endl;

	if ( this->PointFile ) this->ReadIPoints();
	else {
		Timer->StartTimer();
		//Sample : a bit of subsampling (response layer & interpolation)

		//													octave	interval	sample
		FastHessian hessian(this->Integral, this->points, 	4, 		4, 			2, 			this->Threshold);

		if (Mask != 0) hessian.setMask(Mask);
		hessian.getIpoints();
		Timer->StopTimer();
		cout << "FastHessian computed in " << Timer->GetElapsedTime() << "s" << endl;
		/*	Timer->StartTimer();
			hessian.WriteResponseMap();
			Timer->StopTimer();
			cout << "Hessian written in " <<
			Timer->GetElapsedTime() << "s" << endl;
			*/
	}

	Timer->StartTimer();

	if ( this->NumberOfPoints > 0) {

		if ( this->points.size() > this->NumberOfPoints ) {

			partial_sort(this->points.begin(), this->points.begin()+this->NumberOfPoints,
				this->points.end(), compareResponses);
			this->points.resize(this->NumberOfPoints);

		} else {

			sort(this->points.begin(), this->points.end(), compareResponses);

		}

	}

	cout << "Number of keypoints : " << points.size() << endl;

	if ( this->DescriptorType == 0 ) {

		Surf Descriptors( this->Integral, this->points);
		Descriptors.getDescriptors( this->SubVolumeRadius, this->Normalize );

	} else if ( this->DescriptorType == 1 ) {

		Surf Descriptors( this->Integral, this->points);
		Descriptors.getRawDescriptors( this->SubVolumeRadius );

	} else if ( this->DescriptorType == 2 ) {

		vtkNew<vtkMultiThreader> Threader;
		if (this->NbThread > 0) {
			Threader->SetNumberOfThreads(NbThread);
		}
        Threader->SetSingleMethod (ThreadedSubVolumes, (void *) this);
        Timer->StartTimer();
        Threader->SingleMethodExecute ();
    }

	Timer->StopTimer();
	cout << "Descriptors computed in " <<
	Timer->GetElapsedTime() << "s" << endl;
	float min, max;
	min = max = points[ 0 ].descriptor[ 0 ];
	for ( const auto &point : this->points )
		for ( const auto &value : point.descriptor ) {
			if ( max < value ) max = value;
			if ( min > value ) min = value;
		}

	std::cout << "Descriptor range : [" << min << ", " << max << "]" << std::endl;

#ifdef DEBUG
	int rs = 0;
	for(int x=20; x<40; x++)
		for(int y=30; y<60; y++)
			for(int z=40; z<80; z++)
				rs += *static_cast<int*>(this->Cast->GetScalarPointer(x, y, z));

		cout << "rs = " << rs << endl;
		cout << "int= " << BoxIntegral(this->Integral, 20, 30, 40, 20, 30, 40) << endl;
#endif

}


VTK_THREAD_RETURN_TYPE vtk3DSURF::ThreadedSubVolumes (void *arg) {

	vtkMultiThreader::ThreadInfo *Info = (vtkMultiThreader::ThreadInfo*) arg;
	vtk3DSURF *self = (vtk3DSURF *) Info->UserData;
	int MyId = Info->ThreadID;
	int NumberOfThreads = Info->NumberOfThreads;
	const int radius = self->SubVolumeRadius;
	const int size = 8 * radius * radius * radius;
	vtkNew<vtkImageResize> resize;
	resize->SetInputData(self->Resized);
	resize->SetResizeMethodToOutputDimensions ();
	resize->SetOutputDimensions( 2 * radius, 2 * radius, 2 * radius );
	resize->SetNumberOfThreads( 1 );
	resize->SetCropping( 1 );
	double spacing[3];
	self->Resized->GetSpacing( spacing );
	double origin[3];
	self->Resized->GetOrigin( origin );

	for (int id = MyId; id < self->points.size(); id += NumberOfThreads) {

		Ipoint *ipt = &self->points[id];
		float x = ipt->x;
		float y = ipt->y;
		float z = ipt->z;
		float r = radius * ipt->scale;
		float xMin = origin[ 0 ] +  spacing[ 0 ] * ( x - r );
		float xMax = origin[ 0 ] +  spacing[ 0 ] * ( x + r );
		float yMin = origin[ 1 ] +  spacing[ 1 ] * ( y - r );
		float yMax = origin[ 1 ] +  spacing[ 1 ] * ( y + r );
		float zMin = origin[ 2 ] +  spacing[ 2 ] * ( z - r );
		float zMax = origin[ 2 ] +  spacing[ 2 ] * ( z + r );
		resize->SetCroppingRegion( xMin, xMax, yMin, yMax, zMin, zMax );
		resize->Update();
		ipt->allocate( size );
		vtkImageData *output = resize->GetOutput();
		int *voxels = ( int* )  output->GetScalarPointer();
		for ( int i = 0; i < size; i++ ) ipt->descriptor[ i ] = voxels[ i ];

#ifdef DEBUG
		if ( ( ( id % 1000) == 0 ) ) {

			cout << "Point " << id << " done " <<endl;
			cout<< "scale : " << ipt->scale <<endl;;
			cout<< "radius : " << r <<endl;;
			cout<< "ipoint : " <<ipt->x << " "<< ipt->y<< " " << ipt->z<<endl;
			double ori[3];
			double sp[3];
			resize->GetOutput()->GetSpacing(sp );
			cout<<"spaacing : " << sp[0] << " " << sp[1] << " " << sp[2]<< endl;
			resize->GetOutput()->GetOrigin(ori );
			cout<<"Origin: " << ori[0] << " " << ori[1] << " " << ori[2]<< endl;
			vtkNew<vtkMetaImageWriter> writer;
			std::stringstream strfile;
			float coords[3];
			for (int i=0; i < 3; i++) coords[i] =  ori[i] + radius *sp[i];
			cout<< xMin<< " " << xMax << " " << yMin << " " <<yMax << " " << zMin << " " << zMax << endl;
			cout<< "Image center     : " <<coords[0] << " " <<coords[1]<< " " << coords[2]<<endl;
			cout<< "Point coordinates: " << origin[ 0 ] + ipt->x * spacing[ 0 ]
				<< " " << origin[ 1 ] + ipt->y * spacing[ 1 ]
				<< " " << origin[ 2 ] + ipt->z * spacing[ 2 ]<<endl;
			strfile << "Point" << id << ".mhd";
			writer->SetInputData( output );
			writer->SetFileName( strfile.str().c_str() );
			writer->Write();
		}
#endif

	}

	return (VTK_THREAD_RETURN_VALUE);

}


void vtk3DSURF::WritePoints(const char *fileName) {

	double origin[3];
	double spacing[3];
	int dimensions[3];
	this->Cast->GetSpacing(spacing);
	this->Cast->GetOrigin(origin);
	this->Cast->GetDimensions(dimensions);

	picojson::array datapoint;
	for (int i = 0; i != this->points.size(); i++) {
			object iPoint;
			Ipoint &point = this->points[i];
			iPoint["x"] = value(point.x * spacing[ 0 ] + origin[ 0 ]);
			iPoint["y"] = value(point.y * spacing[ 1 ] + origin[ 1 ]);
			iPoint["z"] = value(point.z * spacing[ 2 ] + origin[ 2 ]);
			iPoint["scale"] = value(point.scale * Spacing);
			iPoint["response"] = value (point.response);
			iPoint["laplacian"] = value ((double)point.laplacian);
			picojson::array descriptor;
			bool test = false;
			for (const auto &d : point.descriptor )
				descriptor.push_back( value( d ) );

			iPoint["descriptor"] = value (descriptor);

			datapoint.push_back(value(iPoint));
		}

	picojson::array jorigin;
	picojson::array jspacing;
	picojson::array jdimensions;
	for (int i = 0; i != 3; i++) {
		jorigin.push_back(value(origin[i]));
		jspacing.push_back(value(spacing[i]));
		jdimensions.push_back(value((double)dimensions[i]));
	}

	object root;
	root["origin"] = value(jorigin);
	root["spacing"] = value(jspacing);
	root["dimensions"] = value(jdimensions);
	root["points"] = value(datapoint);

	ofstream pointsFile;
	pointsFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
	pointsFile << value(root);
	pointsFile.close();
}


void vtk3DSURF::WritePointsCSV(const char *fileName) {

	double origin[3];
	double spacing[3];
	int dimensions[3];

	this->Cast->GetSpacing(spacing);
	this->Cast->GetOrigin(origin);
	this->Cast->GetDimensions(dimensions);

	double Sspacing = pow(spacing[0]*spacing[1]*spacing[2], 1.0/3);

	ofstream pointsFile;
	pointsFile.open(fileName, std::ofstream::out | std::ofstream::trunc);
		
	for (int i = 0; i != this->points.size(); i++) {
		Ipoint &point = this->points[i];
		
		pointsFile <<  point.x * spacing[0] + origin[0] << ",";
		pointsFile <<  point.y * spacing[1] + origin[1] << ",";
		pointsFile <<  point.z * spacing[2] + origin[2] << ",";
		pointsFile <<  point.scale * Sspacing           << ",";
		pointsFile <<  point.laplacian                  << ",";
		pointsFile <<  point.response                   << ",";
		for (int k = 0; k < point.descriptor.size(); k++) {
			pointsFile <<  point.descriptor[k];
			if (k<point.descriptor.size() - 1) pointsFile << ",";
		}

		pointsFile << endl;
	}
	pointsFile.close();
}

void vtk3DSURF::WritePointsCSVGZ(const char *fileName, const char *gzOpts) {

	double origin[ 3 ];
	double spacing[ 3 ];
	int dimensions[ 3 ];

	this->Cast->GetSpacing( spacing );
	this->Cast->GetOrigin( origin );
	this->Cast->GetDimensions( dimensions );

	double Sspacing = pow( spacing[ 0 ] * spacing[ 1 ] * spacing[ 2 ], 1.0 / 3.0 );
	string opts( "w" );
	if ( gzOpts ) opts += gzOpts;

	gzFile gz = gzopen( fileName, opts.c_str() );

	for ( int i = 0; i != this->points.size(); i++) {

		Ipoint &point = this->points[ i ];
		gzprintf( gz, "%f,", point.x * spacing[ 0 ] + origin[ 0 ] );
		gzprintf( gz, "%f,", point.y * spacing[ 1 ] + origin[ 1 ] );
		gzprintf( gz, "%f,", point.z * spacing[ 2 ] + origin[ 2 ] );
		gzprintf( gz, "%f,", point.scale * Sspacing );
		gzprintf( gz, "%d,", point.laplacian );
		gzprintf( gz, "%f,", point.response );

		for (int k = 0; k < point.descriptor.size(); k++) {

			if ( k < point.descriptor.size() - 1 ) {

				gzprintf( gz, "%f,", point.descriptor[ k ] );

			} else {

				gzprintf( gz, "%f", point.descriptor[ k ] );

			}

		}

		gzprintf( gz, "\n" );

	}

	gzclose( gz );

}


void vtk3DSURF::WritePointsBinary(const char *fileName) {

	double origin[3];
	double spacing[3];
	int dimensions[3];

	this->Cast->GetSpacing(spacing);
	this->Cast->GetOrigin(origin);
	this->Cast->GetDimensions(dimensions);

	double Sspacing = pow(spacing[0]*spacing[1]*spacing[2], 1.0/3);
	FILE * file = fopen(fileName,"wb");
	float valF;

	for (int i = 0; i != this->points.size(); i++) {
		Ipoint &point = this->points[i];
		valF =  point.x * spacing[0] + origin[0];
		fwrite(&valF, sizeof(float), 1, file);
		valF =  point.y * spacing[1] + origin[1];
		fwrite(&valF, sizeof(float), 1, file);
		valF =  point.z * spacing[2] + origin[2];
		fwrite(&valF, sizeof(float), 1, file);
		valF =  point.scale * Sspacing;
		fwrite(&valF, sizeof(float), 1, file);
		valF = point.laplacian;
		fwrite(&valF, sizeof(float), 1, file);
		valF =  point.response;
		fwrite(&valF, sizeof(float), 1, file);
		fwrite(&point.descriptor[0], sizeof(float), point.descriptor.size(), file);
	}

	fclose(file);
}
