/************************************************************************

   File: InteractorStyleFibreTract.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <vtkRenderWindowInteractor.h>
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include <vtkCommand.h>
#include <vtkWindowLevelLookupTable.h>
#include <vtkActor.h>
#include <vtkMapper.h>

#include "InteractorStyleFibreTract.h"
//#include "otherlibs_include.h"

InteractorStyleFibreTract::InteractorStyleFibreTract()
{
  this->MotionFactor   = 10.0;
  this->_number_of_actors=0;
  _fibre_tract_actor_array=(vtkActor **) malloc(_number_of_actors*sizeof(vtkActor *));
  
}


void InteractorStyleFibreTract::OnRightButtonDown()
{

  if (this->Interactor->GetShiftKey()) 
    {
      _x=this->Interactor->GetEventPosition()[0];
      _y=this->Interactor->GetEventPosition()[1];
      this->StartTractAlphaInteraction();
      
    }
  else
    {
      
      this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
			      this->Interactor->GetEventPosition()[1]);
      if (this->CurrentRenderer == NULL)
	{
	  return;
	}
      
      this->StartDolly();
    }
  


}



void InteractorStyleFibreTract::OnRightButtonUp()
{
  switch (this->State) 
    {
    case VTKIS_DOLLY:
      this->EndDolly();
      break;
    case VTKIS_TRACTALPHA:
      this->EndTractInteraction();
      break;
      
    }
}

void InteractorStyleFibreTract::OnLeftButtonDown()
{

  if (this->Interactor->GetShiftKey()) 
    {
      _x=this->Interactor->GetEventPosition()[0];
      _y=this->Interactor->GetEventPosition()[1];
      this->StartTractColourInteraction();
      
    }
  else
    {
      
      this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], 
			      this->Interactor->GetEventPosition()[1]);
      if (this->CurrentRenderer == NULL)
	{
	  return;
	}
      
      this->StartRotate();
    }
  


}



void InteractorStyleFibreTract::OnLeftButtonUp()
{
  switch (this->State) 
    {
    case VTKIS_ROTATE:
      this->EndRotate();
      break;
    case VTKIS_TRACTCOLOUR:
      this->EndTractInteraction();
      break;
      
    }
}



void InteractorStyleFibreTract::OnMouseMove() 
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
      
    case VTKIS_TRACTALPHA:
      this->ChangeTractAlpha(x-_x,y-_y);
      _x=x;
      _y=y;      
      break;

    case VTKIS_TRACTCOLOUR:
      this->ChangeColourRange(x,y);
      break;

    }
  

}


void InteractorStyleFibreTract::StartTractAlphaInteraction()
{
  
  this->State=VTKIS_TRACTALPHA;
  
}

void InteractorStyleFibreTract::StartTractColourInteraction()
{
  
  this->State=VTKIS_TRACTCOLOUR;
  
}


void InteractorStyleFibreTract::EndTractInteraction()
{
  this->State=VTKIS_NONE;
}

void InteractorStyleFibreTract::AddColourBar(vtkScalarBarActor *colourbar)
{
  m_colourbar=colourbar;
  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor2D(colourbar);

}


//void InteractorStyleFibreTract::ChangeColourRange(float delta_x,float delta_y)
void InteractorStyleFibreTract::ChangeColourRange(float x,float y)
{
  double *range=(double *) malloc(2*sizeof(double));

  int i;
  


  vtkWindowLevelLookupTable *lut;


  for (i=0;i<_number_of_actors;i++)
    {
      


      lut=(vtkWindowLevelLookupTable *) _fibre_tract_actor_array[i]->GetMapper()->GetLookupTable();
  
      //range=lut->GetRange();


      range[0]=x/256;
      range[1]=y/256;

      //range[0]=range[0]+0.1*delta_y;
      //range[1]=range[1]+0.1*delta_x;

      if (range[0]<0.0)
	{
	  range[0]=0.0;
	}
      if (range[1]>1.0)
	{
	  range[1]=1.0;
	}
      if (range[0]>range[1])
	{
	  range[0]=range[1]-0.0001;
	}

      lut->SetRange(range[0],range[1]);

      lut->Build();
  
      cout << range[0] << " " << range[1] << endl;

       	  
      _fibre_tract_actor_array[i]->GetMapper()->SetLookupTable(lut);
      m_colourbar->SetLookupTable(lut);
      m_colourbar->SetNumberOfLabels(10);
    }

  //this changes the colourbar, but it doesn't stay.  The tracts themselves don't change
  this->Interactor->Render();
  //this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->Render();
}


void InteractorStyleFibreTract::ChangeTractAlpha(float delta_x,float delta_y)
{

  double max_rgba[4];
  double min_rgba[4];
  float level;
  float window;
  int i;
  

  vtkWindowLevelLookupTable *lut;


  for (i=0;i<_number_of_actors;i++)
    {
      


      lut=(vtkWindowLevelLookupTable *) _fibre_tract_actor_array[i]->GetMapper()->GetLookupTable();
  

      //4.2:
      //lut->GetMaximumTableValue((double)max_rgba);
      //lut->GetMinimumTableValue((double)min_rgba);

      //4.4 & 5.X
      lut->GetMaximumTableValue(max_rgba[0],max_rgba[1],max_rgba[2],max_rgba[3] );
      lut->GetMinimumTableValue(min_rgba[0],min_rgba[1],min_rgba[2],min_rgba[3] );

  
 

      level=(max_rgba[3]+min_rgba[3])/2;
      
      window=max_rgba[3]-min_rgba[3];
      
      
      level=level+0.01*delta_y;
      //window=window+0.01*delta_x;
      
      
 
      max_rgba[3]=level+window/2;
      min_rgba[3]=level-window/2;
      


   
   
      if (max_rgba[3]<0)
	{
	  max_rgba[3]=0;
	  
	}
      else if (max_rgba[3]>1)
	{
	  max_rgba[3]=1;
	  
	}
      if (min_rgba[3]<0)
	{
	  min_rgba[3]=0;
	  
	}
      else if (min_rgba[3]>1)
	{
	  min_rgba[3]=1;
	  
	}

      lut->SetMaximumTableValue(max_rgba[0],max_rgba[1],max_rgba[2],max_rgba[3]);//or max_rgba
      lut->SetMinimumTableValue(min_rgba[0],min_rgba[1],min_rgba[2],min_rgba[3]);

      lut->Build();
  
  	  
      _fibre_tract_actor_array[i]->GetMapper()->SetLookupTable(lut);
      
    }
  
  this->Interactor->Render();
      
  

}

void InteractorStyleFibreTract::AddActor(vtkActor *fibre_tract_actor)
{
  _number_of_actors+=1;
  
  _fibre_tract_actor_array=(vtkActor **) realloc(_fibre_tract_actor_array,_number_of_actors*sizeof(vtkActor *));
  
  _fibre_tract_actor_array[_number_of_actors-1]=fibre_tract_actor;
  
  
}
