if( CMAKE_PARENT_LIST_FILE STREQUAL CMAKE_CURRENT_LIST_FILE )
  cmake_minimum_required(VERSION 3.30)
  project(vtkOpenSURF3D C CXX)
  find_package(VTK REQUIRED COMPONENTS CommonSystem ImagingColor IOExport )
  include_directories(${VTK_INCLUDE_DIRS})
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED COMPONENTS core)
find_package(ZLIB REQUIRED)
find_package(OpenMP)

option(EXTRACT_CORNER "Extract corner points" OFF)
if(EXTRACT_CORNER)
  add_definitions(-DEXTRACT_CORNER)
endif()

option(BUILD_SURFMATCH "Build surf match program. Requires the TooN library." OFF)

include_directories(${ZLIB_INCLUDE_DIRS})

add_executable(surf3d
  surf3d.cxx
  vtk3DSURF.cxx
  integral.cxx
  fasthessian.cxx
  surf.cxx
)

target_link_libraries(surf3d PRIVATE ${OpenCV_LIBS} ${VTK_LIBRARIES} ${ZLIB_LIBRARIES})
if(TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(surf3d PRIVATE OpenMP::OpenMP_CXX)
elseif(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(BUILD_SURFMATCH)
  add_executable(match3d
    MainMatch.cxx
    MatchPoint.cxx
  )
  target_link_libraries(match3d PRIVATE lapack blas ${VTK_LIBRARIES})
  if(TARGET OpenMP::OpenMP_CXX)
    target_link_libraries(match3d PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
