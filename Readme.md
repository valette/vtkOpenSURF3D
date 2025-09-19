[![CI](https://github.com/valette/vtkOpenSURF3D/actions/workflows/ci.yml/badge.svg)](https://github.com/valette/vtkOpenSURF3D/actions/workflows/ci.yml)

vtkOpenSURF3D : a 3D implementation of Speeded Up Robust Features for VTK
=========================================================================

### Info ###
This code is the SURF3D implementation deriving from this paper:

[1] Rémi Agier, Sébastien Valette, Laurent Fanton, Pierre Croisille and Rémy Prost, Hubless 3D Medical Image Bundle Registration, In Proceedings of the 11th Joint Conference on Computer Vision, Imaging and Computer Graphics Theory and Applications - Volume 3: VISAPP, 265-272, 2016, Rome, Italy.

Which is a 3D extension of the SURF Paper :
[2] Herbert Bay, Andreas Ess, Tinne Tuytelaars, Luc Van Gool, "SURF: Speeded Up Robust Features", Computer Vision and Image Understanding (CVIU), Vol. 110, No. 3, pp. 346--359, 2008

Authors:
* Rémi Agier : most of the code
* Sébastien Valette : minor features, cleanups and packaging
* Chris Evans : writer of the inspiring OpenSURF 2D implementation : https://web.archive.org/web/20150206140535/http://www.chrisevansdev.com/

### Licence ###

GNU GPL

###  Dependencies ###

* CMAKE www.cmake.org
* OpenCV www.opencv.org
* VTK www.vtk.org

###  Compilation Guide ###
	git clone https://github.com/valette/vtkOpenSURF3D.git
	cd vtkOpenSURF3D
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

###  Usage ###

	surf3d file [options]

"file" is a 3D image file. Supported formats are: NIFTI, mhd

Available options:
 * -bin 0/1 : write points as bin file. Default : 0
 * -cmin value : clamp values lower than specified value
 * -cmax value : clamp values larger than specified value
 * -csv 0/1 : write points as csv file. Default : 0
 * -csvgz 0/1 : write points as csv.gz file. Default : 1
 * -json 0/1 : write points as json file. Default : 0
 * -n number : maximum number of points
 * -normalize 0/1 : normalize descriptors (default : 1 )
 * -o basename : set output file name. Default: "points"
 * -r radius : descriptor volume radius. Default : 5
 * -s spacing : resample input image to isotropic sampling with given spacing
 * -t threshold : set detector threshold. Default: 0
 * -type 0/1/2 : descriptor type :
	* 0 : SURF3D descriptor (default). Descriptor size : 48
	* 1 : subvolume HAAR coefficients. Descriptor size : 24 * radius^3
	* 2 : subvolume raw voxels. Descriptor size : 8 * radius^3

### Output keypoints format ###

Each keypoint is described by a list of values:
 * x, y, z coordinates
 * scale
 * sign of laplacian (0 or 1)
 * detector response
 * descriptor (48 components for the SURF3D descriptor)

Note that there is no orientation matrix as we do not estimate orientation (upright SURF) more details here:
https://www.archives-ouvertes.fr/hal-01284240/document

comments, suggestions : https://github.com/valette/vtkOpenSURF3D/issues
