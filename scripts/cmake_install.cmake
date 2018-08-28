# Install script for directory: /home/bic/ilana/links/cvs/mincdiffusion/scripts

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ".")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/ascii_binary_vtk"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/concat-diff-header.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/diff_preprocess.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/dti.m"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/dti_octave.m"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/duplicate_frame"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/extract-diff-frames.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/flip_vector.m"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/grad2FSL.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/header2FSL.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/bvecs_bvals2mnc.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/minc_modify_header-diff.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/mincbet"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/mincresample_vector.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/minctensor.pl"
    "/home/bic/ilana/links/cvs/mincdiffusion/scripts/rgbvector"
    )
endif()

