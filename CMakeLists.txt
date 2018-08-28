PROJECT (MINCDIFFUSION)
cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy) #prevent warnings in cmake 2.6
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#MINC tools
#SET(MINC_BASE   "/usr/local/bic" CACHE PATH "MINC Base dir" )
#SET(NETCDF_BASE "/usr"           CACHE PATH "NetCDF Base dir")
#SET(EZMINC_BASE "ezminc" CACHE PATH "EZMINC Base dir" )
SET(GSL_BASE    "/usr" CACHE PATH  "GSL Base dir" )
#SET(X11_BASE /usr/X11R6 CACHE PATH  "X11 Base dir")
find_package(LIBMINC REQUIRED)
FIND_PACKAGE( X11  REQUIRED )
FIND_PACKAGE( OpenGL REQUIRED )

OPTION(USE_MINC2  "Use MINC2" ON)


#VTK 
FIND_PACKAGE(VTK)
IF (VTK_FOUND)
  INCLUDE (${VTK_USE_FILE})
ENDIF (VTK_FOUND)

SET(CMAKE_CXX_FLAGS -Wno-deprecated)
SET(CMAKE_C_FLAGS   -Wno-deprecated)

SUBDIRS( ezminc/minc_io )


SET(MINCDIFFUSION_C 
  src/misc.c++ src/DisplayWindow.c++
  src/ImageData.c++
  src/Directions.c++ 
  src/DiffusionParameters.c++  
  src/ODFPlotter.c++
  src/S2Interpolator.c++ src/DiffusionODFCalculator.c++
  src/ThreeOrthogonalPlanesPlotter.c++ src/FibreTract.c++
  src/FibreTractSet.c++ src/ProcessedDTIData.c++
  src/FibreTractPlotter.c++ src/LevelSetPlotter.c++
  src/FibreTracking.c++ src/InteractorStyleFibreTract.c++
  src/InteractorStyleMovie.c++ src/InteractorStyleLevelSet.c++
  src/SphericalHarmonics.c++ src/Deconvolver.c++)

  
 
INCLUDE_DIRECTORIES(include ${LIBMINC_INCLUDE_DIRS} ${GSL_BASE}/include/gsl ${GSL_BASE}/include ${NETCDF_BASE}/include ezminc/minc_io ${X11_X11_INCLUDE_PATH} ${OPENGL_INCLUDE_DIR} )


LINK_DIRECTORIES(${LIBMINC_LIBRARY_DIRS} ${GSL_BASE}/lib ${EZMINC_BASE}/lib)

LINK_LIBRARIES(vtkCommon vtkGraphics
  vtkFiltering vtkImaging vtkRendering vtkWidgets vtkIO vtkHybrid 
  minc_io
  m  
  gsl gslcblas 
  glut 
  ${OPENGL_glu_LIBRARY} ${OPENGL_gl_LIBRARY} 
  pthread 
  ${X11_X11_LIB} ${X11_Xext_LIB} 
  ${X11_Xmu_LIB} ${X11_Xt_LIB} ${X11_Xi_LIB} ${X11_SM_LIB} ${X11_ICE_LIB} 
   dl)

IF(USE_MINC2)
  ADD_DEFINITIONS( -DMINC2 )
  LINK_LIBRARIES( ${LIBMINC_LIBRARIES} z)
ELSE(USE_MINC2)
  LINK_LIBRARIES(volume_io minc netcdf)
ENDIF(USE_MINC2)


ADD_LIBRARY(mincdiffusion  ${MINCDIFFUSION_C})
LINK_LIBRARIES(mincdiffusion) #looks like all these programs are using this library

ADD_EXECUTABLE(mincQBI src/QBI.c++ )
ADD_EXECUTABLE(minc3Dvis src/3Dvis.c++ )
ADD_EXECUTABLE(transform_diff_dirs src/transform_diff_dirs.c++ )
ADD_EXECUTABLE(mincFibreTrack src/FibreTrack.c++ )
ADD_EXECUTABLE(mincProbDeconvolve src/ProbDeconvolve.c++)
ADD_EXECUTABLE(mincDeconvolve src/Deconvolve.c++ )
#ADD_EXECUTABLE(mkMIP src/mkMIP.c++)

INSTALL(TARGETS  
  mincDeconvolve 
  mincProbDeconvolve
  mincQBI  
  minc3Dvis 
  transform_diff_dirs 
  mincFibreTrack
  #mkMIP
  RUNTIME DESTINATION bin )

add_subdirectory (scripts)

CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

