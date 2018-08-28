/************************************************************************

   File: InteractorStyleMovie.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <vtkWindowToImageFilter.h>
#include <vtkPNMWriter.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include "InteractorStyleMovie.h"
//#include "otherlibs_include.h"

InteractorStyleMovie::InteractorStyleMovie()
{
  this->MotionFactor   = 10.0;
  m_counter=0;


  sprintf(m_tmpdir,"/tmp/tmp%i",(int) time(NULL));

  sprintf(m_command,"mkdir %s\n",m_tmpdir);
  cout << "tmp dir " << m_tmpdir << endl;
  system(m_command);


  m_window_to_image=vtkWindowToImageFilter::New();
  m_writer=vtkPNMWriter::New();

  //m_window_to_image->SetInput(this->GetInteractor()->GetRenderWindow());
  m_framename=(char *) malloc(200*sizeof(char));
  m_baseframename=(char *) malloc(200*sizeof(char)); // added by IL
  
}


void InteractorStyleMovie::SetBaseName(char * basefilename)
{
  m_basefilename=basefilename;
}

void InteractorStyleMovie::OnLeftButtonDown() 
{ 
  this->FindPokedRenderer(this->GetInteractor()->GetEventPosition()[0], 
                          this->GetInteractor()->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
    {
    return;
    }
  
  if (this->GetInteractor()->GetShiftKey()) 
    {
    if (this->GetInteractor()->GetControlKey()) 
      {
      this->StartDolly();
      }
    else//overriding this so it rotates and saves frames 
      {
      
      this->StartFrameSave();
      }
    } 
  else 
    {
    if (this->GetInteractor()->GetControlKey()) 
      {
      this->StartSpin();
      }
    else 
      {
      this->StartRotate();
      }
    }
}


void InteractorStyleMovie::OnLeftButtonUp()
{
  switch (this->State) 
    {
    case VTKIS_DOLLY:
      this->EndDolly();
      break;

    case VTKIS_PAN:
      this->EndPan();
      break;

    case VTKIS_SPIN:
      this->EndSpin();
      break;

    case VTKIS_ROTATE:
      this->EndRotate();
      break;
    
    case VTKIS_MOVIE:
      this->EndFrameSave();
      break;
    }
}


void InteractorStyleMovie::OnMouseMove() 
{ 
  int x = this->GetInteractor()->GetEventPosition()[0];
  int y = this->GetInteractor()->GetEventPosition()[1];

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
      
    case VTKIS_MOVIE:
      this->SaveFrame();
      this->FindPokedRenderer(x, y);
      this->Rotate();
      this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
      break;
    }
  

}


void InteractorStyleMovie::StartFrameSave()
{
  this->State=VTKIS_MOVIE;
  
  
}



void InteractorStyleMovie::EndFrameSave()
{

  bool ffmpeg=true;

  FILE *fp;
  char mpg_filename[200];
  char param_filename[200];
  this->State=VTKIS_NONE;

  //make the movie:

  sprintf(m_command,"for file in %s/*\ndo\nconvert $file $file.ppm\ndone\n",m_tmpdir);
  //has to be sh, not csh
  //sprintf(m_command,"foreach file (%s/*)\nconvert $file $file.ppm\nend\n",m_tmpdir);

  system(m_command);


  sprintf(mpg_filename,"%s.mp4",m_basefilename); /*mpg seems to bug with player so use mp4*/
  

if (!ffmpeg)
    {

  //write the mpeg_encode param file to tmp dir:

  sprintf(param_filename,"%s/param.txt",m_tmpdir);

  fp=fopen(param_filename,"w"); 

// IL modify input filename to be temp m-baseframename
 fprintf(fp,"PATTERN I\n\nPIXEL HALF\n\nIQSCALE 1\n\nPQSCALE 1\n\nBQSCALE 1\n\nRANGE 10\n\nPSEARCH_ALG EXHAUSTIVE\nBSEARCH_ALG EXHAUSTIVE\n\nOUTPUT %s\n\nGOP_SIZE 1\n\nSLICES_PER_FRAME  10\n\nBASE_FILE_FORMAT PPM\n\nINPUT_CONVERT *\n\nINPUT_DIR %s/\n\nINPUT\n%s*.png.ppm [000-%003d]\nEND_INPUT\n\nREFERENCE_FRAME DECODED\n\n",mpg_filename,m_tmpdir,m_baseframename,m_counter-1); 
  
  fclose(fp);

  //mpeg_encode:

  sprintf(m_command,"mpeg_encode %s\n",param_filename);
  system(m_command);
    }
 else
   {

  /* Use the standard ffmpeg module distributed with Debian
   *  instead of mpeg_encode which requires individual installation
   */
  /*  Note: ffmpeg doesn't seem to like the * for pattern matching, use %03d instead*/ 
   sprintf(m_command,"ffmpeg -y -r 10 -f image2 -i %s/movie%%03d.png.ppm %s", m_tmpdir, mpg_filename);
   system(m_command);
		  
   }
  //rm the files from tmpdir (leaving tmpdir for a future movie) (currently that movie will overwrite the first though: implement prompt for movie name if on second, and keep a counter.
  sprintf(m_command,"rm -f %s/*\n",m_tmpdir);
  system(m_command);


  m_counter=0;

}

void InteractorStyleMovie::SaveFrame()
{

  //not needed: rotated already with Rotate; use this if do an automatic rotate and mk movie: 
  //may need a Modified()?
  //this->GetCurrentRenderer()->GetActiveCamera()->Azimuth( 100 );

  m_window_to_image->SetInput(this->GetCurrentRenderer()->GetRenderWindow());
  m_window_to_image->Modified();
  
  //change name of temp framename so that m_basefilename corresponds to path of output file IL
  sprintf(m_baseframename,"movie",(int) time(NULL));
  sprintf(m_framename,"%s/%s%003d.png",m_tmpdir,m_baseframename,m_counter);
  //sprintf(m_framename,"%s/%s%003d.png",m_tmpdir,m_basefilename,m_counter);
  

  cout << "save frame " << m_framename << endl;
  m_writer->SetInput(m_window_to_image->GetOutput());
  m_writer->SetFileName(m_framename);
  m_writer->Write();
  m_counter+=1;
}

InteractorStyleMovie::~InteractorStyleMovie()
{
  m_writer->Delete();


  sprintf(m_command,"rmdir %s\n",m_tmpdir);
  cout << m_command << endl;

  //system(m_command);

  delete(m_command);
  delete(m_framename);

}
