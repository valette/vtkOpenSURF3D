


find_package(OpenCV REQUIRED)
find_package(VTK REQUIRED)
find_package( ZLIB REQUIRED )

include ("${VTK_USE_FILE}")


SET(NBTHREAD "24" CACHE STRING "number of thread")

option(EXTRACT_CORNER "Extract corner points")

if(EXTRACT_CORNER)
add_definitions( -DEXTRACT_CORNER )
endif()

find_package(OpenMP)

if (OPENMP_FOUND)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)


 
add_definitions( -DNBTHREAD=${NBTHREAD} )

OPTION(BUILD_SURFMATCH "Build surf match program. Requires the Toon library." OFF)

include_directories(${ZLIB_INCLUDE_DIRS})

ADD_EXECUTABLE(surf3d 
 surf3d.cxx  
 vtk3DSURF.cxx
 integral.cxx
 fasthessian.cxx
 ipoint.cxx
 surf.cxx)

#set(CMAKE_CXX_FLAGS "-o3 -pg") 

TARGET_LINK_LIBRARIES(surf3d ${OpenCV_LIBS} ${VTK_LIBRARIES} ${ZLIB_LIBRARIES})

if (BUILD_SURFMATCH)
  ADD_EXECUTABLE(match3d
   MainMatch.cxx
   MatchPoint.cxx
   ipoint.cxx)

  TARGET_LINK_LIBRARIES(match3d lapack blas ${VTK_LIBRARIES})

endif (BUILD_SURFMATCH)
