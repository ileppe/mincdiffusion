/************************************************************************

   File: DisplayWindow.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"

#include "DisplayWindow.h"
#include "otherlibs_include.h"
#include "InteractorStyleFibreTract.h"
#include "InteractorStyleLevelSet.h"
#include "InteractorStyleMovie.h"
#include "LevelSetPlotter.h"



DisplayWindow::DisplayWindow() 
{
  m_window_size=500;
  m_renderwindow = vtkRenderWindow::New();
  m_renderer = vtkRenderer::New();
  m_renderwindowinteractor = vtkRenderWindowInteractor::New();
  m_background_color[0] = 1;
  m_background_color[1] = 1;
  m_background_color[2] = 1;

  //m_renderer->SetBackground(_background_color); 
  m_renderer->SetBackground(m_background_color[0],m_background_color[1], m_background_color[2]); 
  m_renderwindow->AddRenderer(m_renderer);
  m_renderwindow->SetSize(m_window_size,m_window_size);
  m_renderwindowinteractor->SetRenderWindow(m_renderwindow);
  m_renderwindowinteractor->Initialize(); // CLAUDE and ilana  seems you need this not to get abort

  m_interactor_style_level_set_on=false;   //CLAUDE and ilana
  m_interactor_style_fibre_tract_on=false;
  m_interactor_style_movie_on=false;
  m_interactor_style_fibre_tract_on=false;

  //default interactor style is TrackballCamera:  this doesn't work: won't compile - complains about header definition (as if vtkInteractorStyleTrackballCamera
  m_interactor_style_trackball_on=true;
  m_interactor_style_trackball=vtkInteractorStyleTrackballCamera::New();
  m_renderwindowinteractor->SetInteractorStyle(m_interactor_style_trackball);



}


DisplayWindow::~DisplayWindow() 
{

  //DO: Delete all of the Actors! (need an array: need this to delete them while viewing as well)

  //have flags for other (none now) allocated things here.  delete them

  //test for Debian build: DO: figure out how to delete these!!
  //_renderwindow->Delete;
  //_renderer->Delete;
  //_renderwindowinteractor->Delete;

  m_interactor_style_trackball->Delete();


  if (m_interactor_style_level_set_on)
    {      
      m_interactor_style_level_set->Delete();
    }
  else if (m_interactor_style_fibre_tract_on)
    {
      
      m_interactor_style_fibre_tract->Delete();
      
    }
  else if (m_interactor_style_movie_on)
    {
      
      m_interactor_style_movie->Delete();
      
    }
}

void DisplayWindow::AddObject(vtkMapper *Mapper)
{
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(Mapper);
  m_renderer->AddActor(actor);
  m_renderer->Render();
  actor->Delete();
  
}

/*vtkObjext doesn't appear to have GetOutput(): check alt
void DisplayWindow::AddObject(vtkObject *object)
{
  vtkPolyDataMapper *mapper=vtkPolyDataMapper::New();
  mapper->SetInput(object->GetOutput());
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(Mapper);
  m_renderer->AddActor(actor);
  m_renderer->Render();
  actor->Delete();
  mapper->Delete();
}
*/

void DisplayWindow::AddObject(vtkActor *actor)
{
  m_renderer->AddActor(actor);
  //DO: add actor to an internal array of actors (so can delete and re-render)
  m_renderer->Render();
  actor->Delete();
  
}

void DisplayWindow::AddObject2D(vtkActor *actor)
{
  m_renderer->AddActor2D(actor);
  //DO: add actor to an internal array of actors (so can delete and re-render)
  m_renderer->Render();
  actor->Delete();
  
}


void DisplayWindow::AddObject(vtk3DWidget **widget_array, int number_of_widgets)
{
  int i;
 
  for (i=0;i<number_of_widgets;i++) 
    {
     widget_array[i]->SetInteractor(m_renderwindowinteractor);      
     widget_array[i]->On();
    }  
    m_renderer->Render();//IRL
}

void DisplayWindow::AddInteractiveTracts(vtkActor *fibre_tract_actor)
{
  if (!m_interactor_style_fibre_tract_on)
    {
    
      m_interactor_style_fibre_tract=new InteractorStyleFibreTract();
      
      m_renderwindowinteractor->SetInteractorStyle(m_interactor_style_fibre_tract);
      m_interactor_style_fibre_tract_on=true;
      
    }
  
  
  m_interactor_style_fibre_tract->AddActor(fibre_tract_actor);
  
  m_renderer->AddActor(fibre_tract_actor);
  m_renderer->Render();
}

void DisplayWindow::AddInteractiveTracts(vtkActor *fibre_tract_actor, vtkActor *colourbar)
{
  if (!m_interactor_style_fibre_tract_on)
    {
    
      m_interactor_style_fibre_tract=new InteractorStyleFibreTract();
      
      m_renderwindowinteractor->SetInteractorStyle(m_interactor_style_fibre_tract);
      m_interactor_style_fibre_tract_on=true;
      
    }
  
  
  m_interactor_style_fibre_tract->AddActor(fibre_tract_actor);
  m_interactor_style_fibre_tract->AddColourBar((vtkScalarBarActor*) colourbar);
  
  m_renderer->AddActor(fibre_tract_actor);
  m_renderer->Render();
}


void DisplayWindow::MoviePrep(char *moviename)
{
  if (!m_interactor_style_movie_on)
    {
    
      m_interactor_style_movie=new InteractorStyleMovie();
      m_interactor_style_movie->SetBaseName(moviename);
      
      m_renderwindowinteractor->SetInteractorStyle(m_interactor_style_movie);
      m_interactor_style_movie_on=true;
      
    }
  
  m_renderer->Render();
}

void DisplayWindow::AddInteractiveLevelSet(LevelSetPlotter *level_set_plotter)
{
  if (!m_interactor_style_fibre_tract_on)
    {
    
      m_interactor_style_level_set=new InteractorStyleLevelSet();      
      m_renderwindowinteractor->SetInteractorStyle(m_interactor_style_level_set);
      m_interactor_style_level_set_on=true;
      
    }
  
  
  m_interactor_style_level_set->AddPlotter(level_set_plotter,m_renderer);
  
  m_renderer->AddActor(level_set_plotter->GetActor());
  m_renderer->Render();
}

//not implemented
/*
void DisplayWindow::RemoveActor(vtkActor *Actor) 
{
  //not sure if this method would work as is:
  //RemoveActor (method of renderer).  would need to have indices to the actors and store them in this class as some data structure that can add/remove: see Max's visualizeim_with_interface
 
  //_renderer->Render();
}
*/

void DisplayWindow::Display()
{
  m_renderwindowinteractor->Start();
}






/*


keypress user commands should be implemented.  should prompt user where appropriate


will have SetWindowSize, Get, etc., for all variables


  float _rotation[3]={90.0, 90.0, 0.0};
  float _camera_position[3]={0.0, 0.0, 0.0};
  float _camera_focal_point[3]={0.0, 0.0, 0.0};
  float _background_color[3] = {1.0, 1.0, 1.0};


 */

