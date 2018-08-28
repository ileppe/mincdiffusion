/************************************************************************

   File: DisplayWindow.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __DisplayWindow_h
#define __DisplayWindow_h

//#include "otherlibs_include.h"
#include "InteractorStyleFibreTract.h"
#include "InteractorStyleLevelSet.h"
#include "InteractorStyleMovie.h"
#include "LevelSetPlotter.h"
#include "vtkInteractorStyleTrackballCamera.h"

/*
DO: add methods to change initial size of window, etc.  just the defaults (in constructor) for now. 

useful: 'w' converts all actors to wireframe, 's' back to surface, 'e' or 'q' exits, 3 toggles 3D stereo or no, 'p' pust a box around the actor under the mouse, 'r' resets the samera view, 'a' separates the mouse interaction actor by actor: so could move one by one (actor mode), 'c' goes back to moving them all simultaneously  (camera mode), 'f' flys to the point under the cursor, 't' makes the interaction more "grab and drag" trackball, 'j' goes back to joystick, 'u' executes a user-defined function, but I don't know how to define the function yet

could accept arrays of actors and mappers as well.

currently, trackball is the only option.  could remove the default specification of interaction style and then would have vtk's default which has toggling between styles enabled

if map tract point scalars, seg fault on 'q'

currently, windowing of opacity works only for the last tract.  that's because currently, there is only one interctor style and that has only one actor.  needs a re-write: the interactor style needs to have an array of actors and to change all of their luts.

 */

class vtkRenderWindow;
class vtkRenderer;
class vtkRenderWindowInteractor;
class vtkMapper;
class DisplayWindow
{
 public:
  DisplayWindow();
  ~DisplayWindow();
  void AddObject(vtkMapper *Mapper);
  void AddObject(vtkActor *actor);
  void AddObject2D(vtkActor *actor);
  void AddObject(vtk3DWidget **widget_array, int number_of_widgets); 
  //void AddObject(vtkObject *object);
  void AddInteractiveTracts(vtkActor *fibre_tract_actor);
  void AddInteractiveTracts(vtkActor *fibre_tract_actor, vtkActor *colourbar);
  void AddInteractiveLevelSet(LevelSetPlotter *level_set_plotter);
  void Display();
  void CreateMovie(char* filenamebase, int angle);
  void MoviePrep(char *moviename);




 private:  
  vtkInteractorStyleTrackballCamera *m_interactor_style_trackball;
  int m_window_size;
  vtkRenderWindow *m_renderwindow;
  vtkRenderer *m_renderer;
  vtkRenderWindowInteractor *m_renderwindowinteractor;
  double m_background_color[3];
  bool m_interactor_style_fibre_tract_on;
  bool m_interactor_style_movie_on;
  bool m_interactor_style_trackball_on; //not really nec.
  InteractorStyleFibreTract *m_interactor_style_fibre_tract;
  bool m_interactor_style_level_set_on;
  InteractorStyleLevelSet *m_interactor_style_level_set;
  InteractorStyleMovie *m_interactor_style_movie;
  
  
  
};

#endif


