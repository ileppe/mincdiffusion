/************************************************************************

   File: InteractorStyleLevelSet.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkRenderer.h>

#include "InteractorStyleLevelSet.h"
//#include "otherlibs_include.h"
#include "LevelSetPlotter.h"


InteractorStyleLevelSet::InteractorStyleLevelSet()
{
  this->MotionFactor   = 10.0;
  this->_number_of_plotters=0;
  _level_set_plotter_array=(LevelSetPlotter **) malloc(_number_of_plotters*sizeof(LevelSetPlotter *));
  _level_set_renderer_array=(vtkRenderer **) malloc(_number_of_plotters*sizeof(vtkRenderer *));
}


void InteractorStyleLevelSet::OnMiddleButtonDown() 
{
  if (this->Interactor->GetShiftKey()) 
    {
      _x=this->Interactor->GetEventPosition()[0];
      _y=this->Interactor->GetEventPosition()[1];
      this->StartLevelSetInteraction();
      
    }
  else
    {
      
      this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
                          this->Interactor->GetEventPosition()[1]);
      if (this->CurrentRenderer == NULL)
	{
	  return;
	}
  
      this->StartPan();
    }
}



void InteractorStyleLevelSet::OnMiddleButtonUp()
{
  switch (this->State) 
    {
    case VTKIS_PAN:
      this->EndPan();
      break;
    case VTKIS_LEVEL:
      this->EndLevelSetInteraction();
      break;
    }
}



void InteractorStyleLevelSet::OnMouseMove() 
{ 
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  switch (this->State) 
    {
    case VTKIS_ROTATE:
      this->FindPokedRenderer(x, y);
      this->Rotate();
      this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
      break;

    case VTKIS_PAN:
      this->FindPokedRenderer(x, y);
      this->Pan();
      this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
      break;

    case VTKIS_DOLLY:
      this->FindPokedRenderer(x, y);
      this->Dolly();
      this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
      break;

    case VTKIS_SPIN:
      this->FindPokedRenderer(x, y);
      this->Spin();
      this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
      break;
      
    case VTKIS_LEVEL:
      this->ChangeLevel(x-_x,y-_y);
      _x=x;
      _y=y;      
      break;
    }
  

}


void InteractorStyleLevelSet::StartLevelSetInteraction()
{
  this->State=VTKIS_LEVEL;
  
}



void InteractorStyleLevelSet::EndLevelSetInteraction()
{
  this->State=VTKIS_NONE;
}

void InteractorStyleLevelSet::ChangeLevel(float delta_x,float delta_y)
{

  


  int i;
  


  for (i=0;i<_number_of_plotters;i++)
    {
      

      //change high thresh according to delta_x, low according to delta_y  
      //may want to limit low to FLT_EPILSON (because in most cases 0->everything and this might make vis. croak...    
            
      _level_set_plotter_array[i]->SetHighThreshold(_level_set_plotter_array[i]->GetHighThreshold()+0.01*delta_x);
      _level_set_plotter_array[i]->SetLowThreshold(_level_set_plotter_array[i]->GetLowThreshold()+0.01*delta_y);
      
      _level_set_renderer_array[i]->RemoveActor(_level_set_plotter_array[i]->GetActor());
      
      _level_set_renderer_array[i]->AddActor(_level_set_plotter_array[i]->CreateActor());


      
      
    }
  
  this->Interactor->Render();
      
  

}

void InteractorStyleLevelSet::AddPlotter(LevelSetPlotter *level_set_plotter, vtkRenderer *renderer)
{
  _number_of_plotters+=1;
  
  _level_set_plotter_array=(LevelSetPlotter **) realloc(_level_set_plotter_array,_number_of_plotters*sizeof(LevelSetPlotter *));
  
_level_set_renderer_array=(vtkRenderer **) realloc(_level_set_renderer_array,_number_of_plotters*sizeof(vtkRenderer *));

  _level_set_plotter_array[_number_of_plotters-1]=level_set_plotter;
  _level_set_renderer_array[_number_of_plotters-1]=renderer;
  
  
}
