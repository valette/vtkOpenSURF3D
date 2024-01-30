
#include "integral.h"

#include <vtkNew.h>
#include <vtkImageData.h>
#include <vtkTimerLog.h>
//#define DEBUG

void ThreadedIntegral ( int, vtkImageData*, vtkImageData* );

vtkImageData *ComputeIntegral(vtkImageData *source) {

	vtkImageData *integral = vtkImageData::New();

	integral->CopyStructure(source);

	integral->AllocateScalars(VTK_UNSIGNED_LONG_LONG, 1);

	int* dims = integral->GetDimensions();

	ulli* pin_tegral = static_cast<ulli*>(integral->GetScalarPointer());
	int* pso_urce = static_cast<int*>(source->GetScalarPointer());

	long long int* int_egral_incs = integral->GetIncrements();
	long long int* sou_rce_incs = source->GetIncrements();


	vtkNew<vtkTimerLog> Timer;

	Timer->StartTimer();
	ThreadedIntegral( 0, source, integral );
	Timer->StopTimer();
#ifdef DEBUG
	std::cout << "Integral along X done in " << Timer->GetElapsedTime() << "s" << std::endl;
#endif

	Timer->StartTimer();
	ThreadedIntegral( 1, source, integral );
	Timer->StopTimer();
#ifdef DEBUG
	std::cout << "Integral along Y done in " << Timer->GetElapsedTime() << "s" << std::endl;
#endif

	Timer->StartTimer();
	ThreadedIntegral( 2, source, integral );
	Timer->StopTimer();
#ifdef DEBUG
	std::cout << "Integral along Z done in " << Timer->GetElapsedTime() << "s" << std::endl;
#endif


#ifdef DEBUG
	ulli rs = 0;

	for(int x=0; x<dims[0]; x++)
		for(int y=0; y<dims[1]; y++)
			for(int z=0; z<dims[2]; z++)
			{
				rs += sourc(x,y,z);
			}
	cout << "rs  = " << rs << endl;
	cout << "int = " << integ(dims[0]-1,dims[1]-1,dims[2]-1) << endl;

#endif

	return integral;

}


void ThreadedIntegral ( int type, vtkImageData* source, vtkImageData* integral)
{
	ulli* pin_tegral = static_cast<ulli*>(integral->GetScalarPointer());
	int* pso_urce = static_cast<int*>(source->GetScalarPointer());

	long long int* int_egral_incs = integral->GetIncrements();
	long long int* sou_rce_incs = source->GetIncrements();

	int dimensions[3];

	source->GetDimensions(dimensions);

	int dimX = dimensions[0];
	int dimY = dimensions[1];
	int dimZ = dimensions[2];

	if ( type == 2 ) {

		// sum over z
		#pragma omp parallel for
		for (int x = 0; x < dimX; x ++) {
			for (int y = 0; y != dimY; y++) {
				for (int z = 1; z != dimZ; z++) {
					integ(x,y,z) += integ(x,y,z-1);
				}
			}
		}
	} else if ( type == 1 ) {

		// sum over y
		#pragma omp parallel for
		for (int x = 0; x < dimX; x ++) {
			for (int z = 0; z != dimZ; z++) {
				for (int y = 1; y != dimY; y++) {
					integ(x,y,z) += integ(x,y-1,z);
				}
			}
		}
	} else {
		// sum over x
		#pragma omp parallel for
		for (int y = 0; y < dimY; y ++) {
			for (int z = 0; z != dimZ; z++) {
				integ(0,y,z) = sourc(0,y,z);
				for (int x = 1; x != dimX; x++) {
					integ(x,y,z) = integ(x-1,y,z) + sourc(x,y,z);
				}
			}
		}
	}
}
