
#include "integral.h"

#include <vtkImageData.h>
#include <vtkTimerLog.h>
//#define DEBUG


typedef struct {
		vtkImageData *source;
		vtkImageData *integral;
		int ThreadedComputationType;
	 } ThreadData;


vtkImageData *Integral(vtkImageData *source)
{
	vtkImageData *integral = vtkImageData::New();

	integral->CopyStructure(source);

	integral->AllocateScalars(VTK_UNSIGNED_LONG_LONG, 1);

	int* dims = integral->GetDimensions();

	ulli* pin_tegral = static_cast<ulli*>(integral->GetScalarPointer());
	int* pso_urce = static_cast<int*>(source->GetScalarPointer());

	long long int* int_egral_incs = integral->GetIncrements();
	long long int* sou_rce_incs = source->GetIncrements();

/*
	// first Line
	ulli rs = 0;
	for(int z=0; z<dims[2]; z++)
	{
		rs += sourc(0,0,z);
		integ(0,0,z) = rs;
	}

	// first Plane
	for(int y=1; y<dims[1]; y++)
	{
		rs = 0;
		for(int z=0; z<dims[2]; z++)
		{
			rs += sourc(0, y, z);
			integ(0, y, z) = rs + integ(0, y-1, z);
		}
	}

	// Volume ...
	ulli rp[dims[2]]; //Integral image of current x

	for(int x=1; x<dims[0]; x++)
	{
		//Compute first Line of current x
		rs = 0;
		for(int z=0; z<dims[2]; z++)
		{
			rs += sourc(x,0,z);
			rp[z] = rs;
			integ(x, 0, z) = integ(x-1, 0, z) + rs;
		}

		for(int y=1; y<dims[1]; y++)
		{
			rs = 0;
			for(int z=0; z<dims[2]; z++)
			{
				rs += sourc(x,y,z);
				rp[z] += rs; // Attention, contient la derniere ligne explorée (ligne en cours pour [z]<z ligne precédente pour les z par encore atteint)
				integ(x,y,z) = rp[z] + integ(x-1,y,z);
			}
		}
	}
*/

	ThreadData data;
	data.source = source;
	data.integral = integral;

	vtkTimerLog *Timer = vtkTimerLog::New();

	vtkMultiThreader *Threader = vtkMultiThreader::New();
	Threader->SetSingleMethod (ThreadedIntegral, (void *) &data);
	Threader->SetNumberOfThreads(NBTHREAD);

	Timer->StartTimer();
	data.ThreadedComputationType = 0;
	Threader->SingleMethodExecute ();
	Timer->StopTimer();
#ifdef DEBUG
		std::cout << "Integral along X done in " << Timer->GetElapsedTime() << "s" << std::endl;
#endif

	Timer->StartTimer();
	data.ThreadedComputationType = 1;
	Threader->SingleMethodExecute ();
	Timer->StopTimer();
#ifdef DEBUG
		std::cout << "Integral along Y done in " << Timer->GetElapsedTime() << "s" << std::endl;
#endif

	Timer->StartTimer();
	data.ThreadedComputationType = 2;
	Threader->SingleMethodExecute ();
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


VTK_THREAD_RETURN_TYPE ThreadedIntegral (void *arg)
{
	vtkMultiThreader::ThreadInfo *Info = (vtkMultiThreader::ThreadInfo*) arg;
	ThreadData *self = (ThreadData *) Info->UserData;
	int MyId = Info->ThreadID;
	int NumberOfThreads = Info->NumberOfThreads;



	ulli* pin_tegral = static_cast<ulli*>(self->integral->GetScalarPointer());
	int* pso_urce = static_cast<int*>(self->source->GetScalarPointer());

	long long int* int_egral_incs = self->integral->GetIncrements();
	long long int* sou_rce_incs = self->source->GetIncrements();

	int dimensions[3];

	self->source->GetDimensions(dimensions);

	int dimX = dimensions[0];
	int dimY = dimensions[1];
	int dimZ = dimensions[2];








		switch (self->ThreadedComputationType)
	{
	case 2:
		// sum over z
		for (int x = MyId; x < dimX; x += NumberOfThreads) {
			for (int y = 0; y != dimY; y++) {
				for (int z = 1; z != dimZ; z++) {
					integ(x,y,z) += integ(x,y,z-1);
				}
			}
		}
		return (VTK_THREAD_RETURN_VALUE);
	case 1:
		// sum over y
		for (int x = MyId; x < dimX; x += NumberOfThreads) {
			for (int z = 0; z != dimZ; z++) {
				for (int y = 1; y != dimY; y++) {
					integ(x,y,z) += integ(x,y-1,z);
				}
			}
		}
		return (VTK_THREAD_RETURN_VALUE);
	case 0:
	default:
		// sum over x
		for (int y = MyId; y < dimY; y += NumberOfThreads) {
			for (int z = 0; z != dimZ; z++) {
				integ(0,y,z) = sourc(0,y,z);
				for (int x = 1; x != dimX; x++) {
					integ(x,y,z) = integ(x-1,y,z) + sourc(x,y,z);
				}
			}
		}
		return (VTK_THREAD_RETURN_VALUE);
	}
}
