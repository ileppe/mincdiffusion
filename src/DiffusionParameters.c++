/************************************************************************

   File: DiffusionParameters.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include "DiffusionParameters.h"
#include "otherlibs_include.h"
#include "Directions.h"

//using namespace std; //add this for sstream

//not initialized to anything
DiffusionParameters::DiffusionParameters(int n)
{
  _directions=new Directions(n);
  _bvalues=(double *) malloc(n*sizeof(double));
  _number_of_DWIs=n;
  _directions_are_set=true;
  _reference_count=1;
  m_const_b_set=false;
}


DiffusionParameters::DiffusionParameters(const char *minc_filename)
{
  
  minc::minc_1_reader rdr;
  rdr.open(minc_filename,false,true);
        
  std::vector<double> bvalues=rdr.att_value_double("acquisition","bvalues");
        
  _bvalues=new double[bvalues.size()];
        
  std::copy(bvalues.begin(),bvalues.end(),_bvalues);
  
  _directions=new Directions(minc_filename);
  
  _directions_are_set=true;
  _reference_count=1;

}



DiffusionParameters::~DiffusionParameters()
{
  if (_directions_are_set) //they always are
    {
      _directions->Release();
    }
  free(_bvalues);
  //DO: check for more

}

Directions *DiffusionParameters::GetDirections()
{
  return _directions;
  
}

void DiffusionParameters::SetDirections(Directions *directions)
{
  //should check that number is correct!

  if (_directions_are_set)
    {
      _directions->Release();
    }

  _directions=directions;
  _directions->Retain();

  _directions_are_set=true;
  
}



double *DiffusionParameters::GetBValues()
{
  return _bvalues;
  
}

float DiffusionParameters::GetBValue(int n)
{
  //do some checking...n<_number_of_DWIs
  return (float) _bvalues[n];
  
}

void DiffusionParameters::SetBValue(int n, float value)
{
  _bvalues[n]=value;
  
}


int DiffusionParameters::GetConstb()
{
  //they must be used with care: assumes all b values the same, and returns last one
  if (!m_const_b_set)
    {
      m_const_b=this->GetBValue(this->GetNumberOfDWIs()-1);
      m_const_b_set=true;
    }
  return m_const_b;
  
}

void DiffusionParameters::SetConstb(int b)
{
  m_const_b=b;
  m_const_b_set=true;
}

void DiffusionParameters::WriteToMincFile(minc::minc_1_writer& wrt)
{
  

  //should use minc to do this: but currently using minc_modify_header utility:

/*	  
  int n;
  char command[100000];
  strstream commandstream(command,sizeof(command)); //ilana change deprecated header
  
  //ostringstream commandstream;
  commandstream << "minc_modify_header -dinsert acquisition:bvalues=";
  for(n=0;n<this->GetNumberOfDWIs()-1;n++)
    {
      commandstream << this->GetBValue(n) << ",";
    }
  commandstream << this->GetBValue(this->GetNumberOfDWIs()-1) << " ";
  commandstream << minc_filename << endl << ends;
  
  //there's got to be a prettier way to do this
  //string command_str = commandstream.str();
  //const char* command = command_str.c_str();
  
  system(command);
  
  _directions->WriteToMincFile(minc_filename);*/

  std::vector<double> bval(_bvalues,_bvalues+_number_of_DWIs);
  wrt.insert("acquisition","bvalues",bval);
  
   _directions->WriteToMincFile(wrt);
}


void DiffusionParameters::Retain()
{
  _reference_count+=1;
  
}

void DiffusionParameters::Release()
{
  if (_reference_count<2)
    {
      delete(this);
      
    }
  else
    {
      _reference_count-=1;  
    }
  

}


//junk:

/*

  if (n!=_directions->GetNumberofDirections())
    {
      _bvalues_exist=false;
    }


bool DiffusionParameters::BValuesExist()
{
  return _bvalues_exist;
}

bool DiffusionParameters::DirectionsExist()
{
  return _directions->DirectionsExist();
}


 */
