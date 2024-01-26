#ifndef __vtk3DSURF_h
#define __vtk3DSURF_h

#include <vector>
#include <vtkImageData.h>
#include "ipoint.h"

using namespace std;


class vtk3DSURF : public vtkObject 
{

public :
	static vtk3DSURF *New();
	vtkTypeMacro(vtk3DSURF, vtkObject);

	void SetInput(vtkImageData *Input) {
		this->Input = Input;
	};
	
	void SetMask(vtkImageData *Mask) {
		this->Mask = Mask;
	};

	vtkImageData *GetTransformedInput() {
		return this->Cast;
	}
	
	vtkImageData *GetIntegral() {
		return this->IntegralData;
	}

	void Update();

	void WritePoints(const char *fileName);

	void WritePointsCSV(const char *fileName);
	void WritePointsCSVGZ(const char *fileName);
	void WritePointsBinary(const char *fileName);

	vtkSetMacro(Normalize, bool)
	vtkGetMacro(Normalize, bool)

	vtkSetMacro(PointFile, char*)
	vtkGetMacro(PointFile, char*)

	vtkSetMacro(Threshold, double)
	vtkGetMacro(Threshold, double)
	
	vtkSetMacro(NbThread, int)
	vtkGetMacro(NbThread, int)
	
	vtkSetMacro(MaxSize, int)
	vtkGetMacro(MaxSize, int)
	
	vtkSetMacro(Spacing, double)
	vtkGetMacro(Spacing, double)
	
	vtkSetMacro(DescriptorType, int)
	vtkGetMacro(DescriptorType, int)

	vtkSetMacro(NumberOfPoints, int)
	vtkGetMacro(NumberOfPoints, int)

	vtkSetMacro(SubVolumeRadius, int)
	vtkGetMacro(SubVolumeRadius, int)

protected :

	void ReadIPoints();

	//! Get the long descriptors i.e. sub volumes
	static VTK_THREAD_RETURN_TYPE ThreadedSubVolumes (void *arg);

	vector<Ipoint> points;

	vtkImageData *Input;

	vtkImageData *Mask;
	
	vtkImageData *Cast;
	vtkImageData *Resized;
	
	vtkImageData *IntegralData;

	double Threshold;
	int NbThread;
	int MaxSize;
	double Spacing;
	int DescriptorType;
	int SubVolumeRadius;
	int NumberOfPoints;
	char *PointFile;
	bool Normalize;

	/// the constructor
	vtk3DSURF() {
		this->Cast = 0;
		this->IntegralData = 0;
		this->Mask = 0;
		this->Normalize = true;
		this->Threshold = 0.0004;
		this->NbThread = -1;
		this->MaxSize = 100;
		this->DescriptorType = 0;
		this->NumberOfPoints = -1;
		this->PointFile = 0;
		this->Resized = 0;
		this->Spacing = 0;
		this->SubVolumeRadius = 5;
	}

	/// the destructor
	 ~vtk3DSURF() {
		if (this->IntegralData) {
			this->IntegralData->Delete();
		}
		if (this->Cast) {
			this->Cast->Delete();
		}

		if (this->Mask) {
			this->Mask->Delete();
		}
	}

};
#endif
