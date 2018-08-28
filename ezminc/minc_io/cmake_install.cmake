# Install script for directory: /home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/libminc_io.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_io_exceptions.h"
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_io_fixed_vector.h"
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_io_simple_volume.h"
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_1_rw.h"
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_1_simple.h"
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_1_simple_rw.h"
    "/home/bic/ilana/links/cvs/mincdiffusion/ezminc/minc_io/minc_io_4d_volume.h"
    )
endif()

