/************************************************************************

   File: FibreTractPlotter.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __FibreTractPlotter_h
#define __FibreTractPlotter_h 

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __FibreTractSet_h
#include "FibreTractSet.h"
#endif



/*

the colour bar for colouring tracts by scalars is wrong for the 0.1ish values: they display as violet!!!

tubes are not implemented for thresholding: not sure about how, and would be too slow if a lot of data

thresholding of indices not done

show tracts above a certain threshold of scalar index in file. 

DO: could have multiple scalar indices in file (currently 1), in which case would choose one here.  or, better still, toggle in vis

try using vtkSmoothPolyDataFilter.h to smooth the tracts (else splines)

Note that this filter operates on the lines, polygons, and triangle strips composing an instance of vtkPolyData. Vertex or poly-vertex cells are never modified. 

this seems to always complain about no normals to line and then plot anyway

the color 0 0 0 is invisible (check): just set the color to this (or do transparency as before) below threshold

DO: enable radius varying by diffusion (or fibre) ODF value tangent to tract at each point: would need to save these with the tracts: a lot of data..

DO; add smoothing of tracts: splines.  not sure to what extent if any this is already done by the mapper (I think it is not)  

 */


class FibreTractPlotter
{
 public:
  FibreTractPlotter(char *tracts_filename);
  ~FibreTractPlotter();
  void SetTubesOn();
  void SetTubesOn(float tube_radius);
  void SetColour(double red, double green, double blue); /*change to double for vtk 5 ilana*/
  void MapPointScalarsOn();
  //void MapPointRGBOn();
  vtkActor *CreateActor();
  void ColourByDirection();
  void ColourByScalarFile(char* tract_scalars_name, float low, float high);
  vtkActor *ColourBar();

 private:
  bool _plot_as_tubes;
  float _tube_radius;
  //FibreTractSet _fibre_tract_set;
  vtkPolyData *_tracts;
  double _red; /*change to double for vtk 5 ilana*/
  double _green;
  double _blue;
  bool _map_point_scalars;
  bool _map_point_rgb;
  //tract_scalars?
  bool m_colour_by_direction;
  //bool m_colour_by_scalars;
  vtkLookupTable *m_lut;
  char *m_tracts_filename;
  bool m_map_file_scalars;
  float m_file_low, m_file_high;
};


#endif
