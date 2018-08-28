/************************************************************************

   File: InteractorStyleFibreTract.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __InteractorStyleFibreTract_h
#define __InteractorStyleFibreTract_h

//#include <vtkInteractorStyle.h>
#include <vtkScalarBarActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
//#include "otherlibs_include.h"
#define VTKIS_TRACTALPHA 10
#define VTKIS_TRACTCOLOUR 11
/*

could consider changing or adding to this so that have a hard threshold on the tracts displayed instead of smooth falloff to fuzzy

could leave those below thresh at at translucent and those above opaque

all tracts in a scene should have the SAME lut.  i.e. the same mapping from number to alpha.  could do this by taking the lut for first, chaning it, and getting the lut fr all the others, getting the color, setting the lut equal to the 1st tract set's lut, then changing the min/max color.
 */

class InteractorStyleFibreTract : public vtkInteractorStyleTrackballCamera
{
 public:
  InteractorStyleFibreTract();
  
  void OnRightButtonDown();
  void OnRightButtonUp();
  void OnLeftButtonDown();
  void OnLeftButtonUp();
  void OnMouseMove();
  void StartTractAlphaInteraction();
  void StartTractColourInteraction();
  void EndTractInteraction();
  void ChangeColourRange(float x,float y);
  void ChangeTractAlpha(float delta_x,float delta_y);
  void AddActor(vtkActor *fibre_tract_actor);
  void AddColourBar(vtkScalarBarActor *colourbar);


 private:
  float _x;
  float _y;
  vtkActor **_fibre_tract_actor_array;
  int _number_of_actors;
  vtkScalarBarActor *m_colourbar;
  
};


#endif
