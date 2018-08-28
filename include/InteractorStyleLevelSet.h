/************************************************************************

   File: InteractorStyleLevelSet.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __InteractorStyleLevelSet_h
#define __InteractorStyleLevelSet_h

#include <vtkInteractorStyleTrackballCamera.h>
//#include "otherlibs_include.h"
#include "LevelSetPlotter.h"


#define VTKIS_LEVEL 10

/*

could have a fuzzy outline or something here

should display the threshold 

could: not let low go to 0 or something..

for speed, could pre-calculate the surfaces and store them and them just flip through

 */
 
class vtkRenderer;
class LevelSetPlotter;
class InteractorStyleLevelSet : public vtkInteractorStyleTrackballCamera
{
 public:
  InteractorStyleLevelSet();
  
  void OnMiddleButtonDown();
  void OnMiddleButtonUp();
  void OnMouseMove();
  void StartLevelSetInteraction();
  void EndLevelSetInteraction();
  void ChangeLevel(float delta_x,float delta_y);
  void AddPlotter(LevelSetPlotter *level_set_plotter, vtkRenderer *renderer);
  


 private:
  float _x;
  float _y;
  LevelSetPlotter **_level_set_plotter_array;
  vtkRenderer **_level_set_renderer_array;
  int _number_of_plotters;
  
  
};


#endif
