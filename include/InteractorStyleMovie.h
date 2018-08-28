/************************************************************************

   File: InteractorStyleMovie.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __InteractorStyleMovie_h
#define __InteractorStyleMovie_h

//#include "otherlibs_include.h"
#include <vtkInteractorStyleTrackballCamera.h>

#define VTKIS_MOVIE 10

//could add an opiton to automically rotate through a user-defined angle too (keystroke initiated, give angle in 3Dvis along with filename).

//could make this more all-encompassing: *anything* + shift-> movie (just rework more functions from TrackBallCamera.

//possible that we want this to be Joystick instead (easy to change)


//may want a -clobber option for the output .mpg, and an option for multiple movies (prompt for name each time...)  (can rename already written one for now)

class vtkWindowToImageFilter;
class vtkPNMWriter;
class InteractorStyleMovie : public vtkInteractorStyleTrackballCamera
{
 public:
  InteractorStyleMovie();
  void SetBaseName(char * basefilename);
  void OnLeftButtonDown();
  void OnLeftButtonUp();
  void OnMouseMove();
  ~InteractorStyleMovie();


 private:
  void SaveFrame();
  void StartFrameSave();
  void EndFrameSave();
  char *m_basefilename;
  int m_counter;
  vtkWindowToImageFilter *m_window_to_image;
  vtkPNMWriter *m_writer;
  char *m_framename, *m_baseframename;
  char m_command[300];
  char m_tmpdir[200];


};


#endif
