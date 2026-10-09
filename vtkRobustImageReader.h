#ifndef __vtkRobustImageReader_h
#define __vtkRobustImageReader_h

#include <iostream>
#include <sstream>

#include <vtkImageData.h>
#include <vtkImageFlip.h>
#include <vtkImageReader2Factory.h>
#include <vtkImageShiftScale.h>
#include <vtkMatrix4x4.h>
#include <vtkMetaImageReader.h>
#include <vtkNew.h>
#include <vtkNIFTIImageReader.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>

// v3 : read both QForm and SForm matrices
// v2 : shift and scale values read by nifti reader when needed
// v1


class vtkRobustImageReader : public vtkObject
{

public :

	static vtkRobustImageReader *New();

	vtkTypeMacro(vtkRobustImageReader,vtkObject);

	void Update() {
		if (!this->FileName) {
			vtkErrorMacro(<< "FileName not set");
			return;
		}

		vtkNew<vtkImageReader2Factory> imageReaderFactory;
		vtkNew<vtkMetaImageReader> metaImageReader;
		imageReaderFactory->RegisterReader(metaImageReader);
		vtkNew<vtkNIFTIImageReader> niftiiImageReader;
		imageReaderFactory->RegisterReader(niftiiImageReader);

		// Create a reader for image and try to load it
		vtkSmartPointer<vtkImageReader2> Reader =
			vtkSmartPointer<vtkImageReader2>::Take(imageReaderFactory->CreateImageReader2(FileName));

		if (!Reader) {
			vtkErrorMacro(<< "Cannot load file " << FileName << " as an image file.");
			this->Output = nullptr;
			return;
		}

		Reader->SetFileName(FileName);
		Reader->Update();
		Output = Reader->GetOutput();

		bool flip[3] = {false, false, false};
		double Origin[3] = {0.0, 0.0, 0.0};
		double Spacing[3];
		int Dimensions[3];

		if (strcmp (Reader->GetClassName(), "vtkMetaImageReader") == 0) {
			std::ifstream is(FileName);
			std::string line, prefix;
			double transformationMatrix[9] ;

			while(is) {
				getline(is, line);
				std::istringstream iss(line);
				iss >> prefix;
				if ((prefix=="TransformMatrix") || (prefix == "Orientation") || (prefix == "Rotation")) {
					iss >> prefix; // eat the "=" sign
					for (int i = 0; i < 9; i++) {
						iss >> transformationMatrix[i];
					}

					for (int i = 0; i < 3; i++) {
						if (transformationMatrix[4 * i] < 0) {
							flip[i] = true;
						}
					}
				}
			}
		}

		if (strcmp (Reader->GetClassName(), "vtkNIFTIImageReader") == 0) {
		    vtkNIFTIImageReader *niftiReader = (vtkNIFTIImageReader*) Reader.GetPointer();

            vtkSmartPointer<vtkMatrix4x4> defaultMatrix;
            vtkMatrix4x4 *formMatrix = niftiReader->GetQFormMatrix();
			if (!formMatrix) {
				formMatrix = niftiReader->GetSFormMatrix();
				if (!formMatrix) {
					defaultMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
					formMatrix = defaultMatrix;
				}
            }

            for (int i = 0; i < 3; i++) {
				if (formMatrix->GetElement(i, i) < 0) {
					flip[i] = true;
				}
				Origin[i] = formMatrix->GetElement(i,3);
			}


            if (niftiReader->GetQFac() < 0) {
				for (int i = 0; i < 2; i++) {
					flip[i] = !flip[i];
					Origin[i] = -Origin[i];
				}
            }

			Output->SetOrigin(Origin);
		}

		Output->GetSpacing(Spacing);
		Output->GetDimensions(Dimensions);

        if (flip[0] || flip[1] || flip[2]) {
            for (int i = 0; i < 3; i++) {
                if (!flip[i]) continue;
                std::cout << "Warning : RobustReader flipping dimension " << i << std::endl;
                vtkNew<vtkImageFlip> imgFlip;
                imgFlip->SetInputData(Output);
                imgFlip->SetFilteredAxis (i);
                imgFlip->Update();
                Output = imgFlip->GetOutput();
				Output->GetOrigin(Origin);
                Origin[i] = Origin[i] - Spacing[i] * ( Dimensions[i] - 1);
				Output->SetOrigin(Origin);
            }
        }

		if (strcmp (Reader->GetClassName(), "vtkNIFTIImageReader") == 0) {
		    vtkNIFTIImageReader *niftiReader = (vtkNIFTIImageReader *) Reader.GetPointer();

			double slope = niftiReader->GetRescaleSlope();
			double intercept = niftiReader->GetRescaleIntercept();

			if ( ( slope != 1.0 ) || ( intercept != 0 ) ) {

				vtkNew<vtkImageShiftScale> shiftScale;
				shiftScale->SetShift( intercept );
				shiftScale->SetScale( slope );
				shiftScale->SetInputData( Output );
				shiftScale->Update();
				Output = shiftScale->GetOutput();
				std::cout << "Warning : RobustReader shifting nifti volume values" << std::endl;

			}

		}

	}

	vtkGetObjectMacro(Output, vtkImageData)

	vtkGetStringMacro(FileName)
	vtkSetStringMacro(FileName)

protected :
	vtkSmartPointer<vtkImageData> Output;

	char *FileName;

	vtkRobustImageReader() {
		Output = nullptr;
		FileName = nullptr;
	}

	~vtkRobustImageReader() override {
		this->SetFileName(nullptr);
	}
};

#endif

vtkStandardNewMacro(vtkRobustImageReader);

