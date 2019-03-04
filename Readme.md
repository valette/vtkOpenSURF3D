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

* VTK www.vtk.org
* CMAKE www.cmake.org

###  Compilation Guide ###
	git clone https://github.com/valette/vtkOpenSURF3D.git
	cd vtkOpenSURF3D
	cmake . -DCMAKE_BUILD_TYPE=Release
	make

comments, suggestions : https://github.com/valette/vtkOpenSURF3D/issues