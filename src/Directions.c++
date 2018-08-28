/************************************************************************

   File: Directions.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <float.h>

#include "Directions.h"
#include "otherlibs_include.h"
//#include <minc.h>
#include "misc.h"

#include <volume_io.h>
//using namespace std; //add this for new sstream

Directions::Directions(int n)
{
  
  int i;
  
  if (n>0)
    {
      for (i=0;i<3;i++) 
	{
	  m_directions[i]=(double *) malloc(n*sizeof(double));
	}
    }
  m_number_of_directions=n;
  m_reference_count=1;
  
}

Directions::Directions(const char *filename)
{
  int i;
  std::string fname(filename);
  
  
  if (fname.rfind(".mnc")==(fname.length()-4))
  {
	
	/*
	int cdfid;
	int acquisition_varid;
	int len;
	
	
	cdfid=ncopen(filename,NC_NOWRITE);
	
	acquisition_varid=ncvarid(cdfid,"acquisition");
	
	//these (ncattinq and ncattget) might return errors: they return int: check:
	ncattinq(cdfid,acquisition_varid,"direction_x",NULL,&len);
	
	//use len to allocate directions: assume double
	
	//currently allocating direction_y and direction_z to be the same size.  DO: get these the same way and check that they are all the same.
	
	for (i=0;i<3;i++) 
	   {
	  	m_directions[i]=(double *) malloc(len*sizeof(double));
	   }
	
	
	ncattget(cdfid,acquisition_varid,"direction_x",m_directions[0]);
	ncattget(cdfid,acquisition_varid,"direction_y",m_directions[1]);
	ncattget(cdfid,acquisition_varid,"direction_z",m_directions[2]);*/
        minc::minc_1_reader rdr;
        rdr.open(filename,false,true);
        
        std::vector<double> dir_x=rdr.att_value_double("acquisition","direction_x");
        std::vector<double> dir_y=rdr.att_value_double("acquisition","direction_y");
        std::vector<double> dir_z=rdr.att_value_double("acquisition","direction_z");
        
        m_directions[0]=new double[dir_x.size()];
        m_directions[1]=new double[dir_y.size()];
        m_directions[2]=new double[dir_z.size()];
        
        std::copy(dir_x.begin(),dir_x.end(),m_directions[0]);
        std::copy(dir_y.begin(),dir_y.end(),m_directions[1]);
        std::copy(dir_z.begin(),dir_z.end(),m_directions[2]);
		
	/* DO: check error handling from above.  when error:
	   {
	   Error("Error reading directions from minc file"); 
	   }
	*/
	
	m_number_of_directions=dir_x.size();
	
	//ncclose(cdfid);
  } else  { //assume text file with directions
	  
	  int numdirs=0;
	  int allocated;
	  int alloc_factor=500;
	  

	  ifstream fd(filename);
	  if (!fd)
	  {
	  	Error("cannot open file");
	  }
	  
	  //allocate alloc_factor directions initially and realloc if necessary:
	  for (i=0;i<3;i++) 
	    {
	      m_directions[i]=(double *) malloc(alloc_factor*sizeof(double));
	    }
	  allocated=alloc_factor;
	  

	  while(fd >> m_directions[0][numdirs] >> m_directions[1][numdirs] >> m_directions[2][numdirs])
	    {		
	      
	      numdirs+=1;
	      //cout << "new dir " << m_directions[0][numdirs-1] << " " << numdirs-1 << endl;
	      if (numdirs>allocated)
		{
		  for (i=0;i<3;i++) 
		    {
		      m_directions[i]=(double *) realloc(m_directions[i],(allocated+alloc_factor)*sizeof(double));
		    } 
		  allocated+=alloc_factor;
		  
		}
	    }
	  
	  m_number_of_directions=numdirs;
	}
  
  m_reference_count=1;
}


Directions::~Directions()
{
  int i;
  
  for (i=0;i<3;i++)
    {
      free(m_directions[i]);
      
    }
  
}


int Directions::GetNumberOfDirections()
{

  return m_number_of_directions;
  
}


float Directions::GetXComponent(int n)
{
  return (float) m_directions[0][n];
  
}

float Directions::GetYComponent(int n)
{
  return (float) m_directions[1][n];
}


float Directions::GetZComponent(int n)
{
  return (float) m_directions[2][n];
}


void Directions::SetXComponent(int n, float xcomp)
{
  m_directions[0][n]=(double) xcomp;
}



void Directions::SetYComponent(int n, float ycomp)
{
  m_directions[1][n]=(double) ycomp;
}


void Directions::SetZComponent(int n, float zcomp)
{
  m_directions[2][n]=(double) zcomp;
}

VECTOR3D Directions::GetDirection(int n)
{
  VECTOR3D vector;
  
  vector.x=(float) m_directions[0][n];
  vector.y=(float) m_directions[1][n];
  vector.z=(float) m_directions[2][n];

  return vector;
  
}

void Directions::SetDirection(int n, VECTOR3D vector)
{
  //DO: if n>number of directions, realloc or return error or warning
  m_directions[0][n]=(double) vector.x;
  m_directions[1][n]=(double) vector.y;
  m_directions[2][n]=(double) vector.z;
}




void Directions::WriteToMincFile(minc::minc_1_base& wrt)
{
  //DO: write to header using ncattput instead of using system call to minc_modify_header: faster and less clunky.  However, I get error "ncvarid: ncid 11: Variable not found" with this option: undiagnosed at this point

/*
  int n;
  //Q:how to I do this dynamically??  also, could forgo the whole commandstream and just strcat
  char command[100000];
  strstream commandstream(command,sizeof(command)); //ilana change deprecated header
  //ostringstream commandstream;

  int acquisition_varid,cdfid;

  bool minc_modify_header=true;


  if (minc_modify_header)
    {
  commandstream << "minc_modify_header -dinsert acquisition:direction_x=";
  for(n=0;n<this->GetNumberOfDirections()-1;n++)
    {
      commandstream << this->GetXComponent(n) << ",";
    }
  commandstream << this->GetXComponent(this->GetNumberOfDirections()-1) << " ";
  commandstream << minc_filename << endl;
  
  commandstream << "minc_modify_header -dinsert acquisition:direction_y=";
  for(n=0;n<this->GetNumberOfDirections()-1;n++)
    {
      commandstream << this->GetYComponent(n) << ",";
    }
  commandstream << this->GetYComponent(this->GetNumberOfDirections()-1) << " ";
  commandstream << minc_filename << endl;
  
  commandstream << "minc_modify_header -dinsert acquisition:direction_z=";
  for(n=0;n<this->GetNumberOfDirections()-1;n++)
    {
      commandstream << this->GetZComponent(n) << ",";
    }
  commandstream << this->GetZComponent(this->GetNumberOfDirections()-1) << " ";
  commandstream << minc_filename << endl << ends;
  
  //there's got to be a prettier way to do this
  //string command_str = commandstream.str();
  //const char* command = command_str.c_str();

  system(command);
    }

  else //use netcdf directly:
    {
      



      cdfid=ncopen(minc_filename,NC_WRITE);

      ncredef(cdfid);
      acquisition_varid=ncvarid(cdfid,"acquisition");
	
      
      ncattput(cdfid,acquisition_varid,"direction_x",NC_FLOAT,m_number_of_directions,m_directions[0]);
      ncattput(cdfid,acquisition_varid,"direction_y",NC_FLOAT,m_number_of_directions,m_directions[1]);
      ncattput(cdfid,acquisition_varid,"direction_z",NC_FLOAT,m_number_of_directions,m_directions[2]);
      ncendef(cdfid); 
      
      ncclose(cdfid);
    }
  
*/
  std::vector<double> dir_x(m_directions[0],m_directions[0]+m_number_of_directions);
  std::vector<double> dir_y(m_directions[1],m_directions[1]+m_number_of_directions);
  std::vector<double> dir_z(m_directions[2],m_directions[2]+m_number_of_directions);
  wrt.insert("acquisition","direction_x",dir_x);
  wrt.insert("acquisition","direction_y",dir_y);
  wrt.insert("acquisition","direction_z",dir_z);
  

}




/*void Directions::WriteToMincFile(char *minc_filename, bool clobber)
{
  bool may_output=true;
  char use_new_filename;
  char new_minc_filename[1000];

  strcpy(new_minc_filename,minc_filename);



  while (may_output)
    {
      if (!clobber)
	{
	  clobber=check_clobber_file(new_minc_filename);
	}
  
      if (clobber)
	{
	  this->WriteToMincFile(new_minc_filename);
	  may_output=false;
	  
	  
	}
      else //prompt user for new filename
	{
	  cout << "Do you wish to use an alternate file name? [y] [n]\n";
	  cin >> use_new_filename;
	  if (use_new_filename=='y')
	    {
	      cout << "new file name:\n";
	      cin >> new_minc_filename;
	    }
	  else
	    {
	      may_output=false;
	      
	    }
	  
	}
    }
  
      
  

}*/

Directions *Directions::Transform(float rotation[3][3])
{
  int i;
  Directions *new_directions=new Directions(m_number_of_directions);

  for (i=0;i<m_number_of_directions;i++) 
    {
      new_directions->SetXComponent(i,rotation[0][0]*this->GetXComponent(i)+rotation[0][1]*this->GetYComponent(i)+rotation[0][2]*this->GetZComponent(i));
      new_directions->SetYComponent(i,rotation[1][0]*this->GetXComponent(i)+rotation[1][1]*this->GetYComponent(i)+rotation[1][2]*this->GetZComponent(i));
      new_directions->SetZComponent(i, rotation[2][0]*this->GetXComponent(i)+rotation[2][1]*this->GetYComponent(i)+rotation[2][2]*this->GetZComponent(i));
      

    }
  
  //if truly a rotation, this won't be nec., but if a scaling too, it is
  new_directions->Normalize();
  
  return new_directions;
}

Directions *Directions::Transform(const char *xfm_filename)
{
  int i;
  
  //read in the transform file
  VIO_General_transform general_transform;

  if(input_transform_file((char*)xfm_filename, &general_transform)!=VIO_OK) 
    Error("Can't open XFM file");


  ::VIO_Transform *transform=get_linear_transform_ptr(&general_transform);

  //check Transform type to get at the components and extract the non-shift part 
  
  //create new directions and set all components equal to transformed (and re-normalized)

  Directions *new_directions=new Directions(m_number_of_directions);


  //get transform values:

  VIO_Real tx[3];
  VIO_Real ty[3];  
  VIO_Real tz[3];

  get_transform_x_axis_real(transform,tx);
  get_transform_y_axis_real(transform,ty);
  get_transform_z_axis_real(transform,tz);
    

  for (i=0;i<m_number_of_directions;i++) 
    {

      new_directions->SetXComponent(i,(float) tx[0]*this->GetXComponent(i)+ty[0]*this->GetYComponent(i)+tz[0]*this->GetZComponent(i));
      new_directions->SetYComponent(i,(float) tx[1]*this->GetXComponent(i)+ty[1]*this->GetYComponent(i)+tz[1]*this->GetZComponent(i));
      new_directions->SetZComponent(i,(float) tx[2]*this->GetXComponent(i)+ty[2]*this->GetYComponent(i)+tz[2]*this->GetZComponent(i));



    }
  
  //renormalize:

  new_directions->Normalize();
  
  return new_directions;
}

void Directions::IncreaseNumberOfDirections(int n)
{
  int i;
  
  if (this->GetNumberOfDirections()>0)
    {
      for (i=0;i<3;i++) 
	{
	  
	  m_directions[i]=(double *) realloc(m_directions[i],(m_number_of_directions+n)*sizeof(double));
	 
	}
    }
  else
    {
      for (i=0;i<3;i++) //note that realloc can be called with NULL pointer and behave as malloc
	{
	  m_directions[i]=(double *) malloc(n*sizeof(double));
	}
    }
  m_number_of_directions+=n;
 

}

vtkActor *Directions::GetActor()
{
  //could do this as lines instead; could add a translucent sphere
  int i;
  
  vtkPolyData *polydata=vtkPolyData::New();
  vtkPoints *points=vtkPoints::New();

  for (i=0;i<m_number_of_directions;i++) 
    {
      points->InsertPoint(i,this->GetXComponent(i),this->GetYComponent(i),this->GetZComponent(i));
    }

  polydata->SetPoints(points);   
  vtkSphereSource *sphere=vtkSphereSource::New();
  float radius=0.02;
  sphere->SetRadius(radius);
  vtkGlyph3D *spheres=vtkGlyph3D::New();
  spheres->SetInput(polydata);
  spheres->SetSource(sphere->GetOutput());
  vtkPolyDataMapper *glyphmapper=vtkPolyDataMapper::New();
  glyphmapper->SetInput(spheres->GetOutput());
  vtkActor *glyphactor=vtkActor::New();
  glyphactor->SetMapper(glyphmapper);
  glyphactor->GetProperty()->SetColor(0,0,1);


  spheres->Delete();
  sphere->Delete();
  polydata->Delete();
  points->Delete();

  return glyphactor;
  
}

vtkActor *Directions::GetActor(float red, float green, float blue)
{
  //could do this as lines instead; could add a translucent sphere
  int i;
  
  vtkPolyData *polydata=vtkPolyData::New();
  vtkPoints *points=vtkPoints::New();

  for (i=0;i<m_number_of_directions;i++) 
    {
      points->InsertPoint(i,this->GetXComponent(i),this->GetYComponent(i),this->GetZComponent(i));
    }

  polydata->SetPoints(points);   
  vtkSphereSource *sphere=vtkSphereSource::New();
  float radius=0.02;
  sphere->SetRadius(radius);
  vtkGlyph3D *spheres=vtkGlyph3D::New();
  spheres->SetInput(polydata);
  spheres->SetSource(sphere->GetOutput());
  vtkPolyDataMapper *glyphmapper=vtkPolyDataMapper::New();
  glyphmapper->SetInput(spheres->GetOutput());
  vtkActor *glyphactor=vtkActor::New();
  glyphactor->SetMapper(glyphmapper);
  glyphactor->GetProperty()->SetColor(red,green,blue);


  spheres->Delete();
  sphere->Delete();
  polydata->Delete();
  points->Delete();

  return glyphactor;
  
}

void Directions::Normalize()
{
  int i;
  float mag;
  
  
  for (i=0;i<m_number_of_directions;i++)
    {
      mag=sqrt(pow(m_directions[0][i],2)+pow(m_directions[1][i],2)+pow(m_directions[2][i],2));
      if (mag>FLT_EPSILON)
	{
	  m_directions[0][i]=m_directions[0][i]/mag;
	  m_directions[1][i]=m_directions[1][i]/mag;
	  m_directions[2][i]=m_directions[2][i]/mag;
	  
	}
      
      
    }
  

}

void Directions::Normalize(int n)
{
  float mag;
  

  mag=sqrt(pow(m_directions[0][n],2)+pow(m_directions[1][n],2)+pow(m_directions[2][n],2));
  if (mag>FLT_EPSILON)
    {
      m_directions[0][n]=m_directions[0][n]/mag;
      m_directions[1][n]=m_directions[1][n]/mag;
      m_directions[2][n]=m_directions[2][n]/mag;
      
    }

}


void Directions::Retain()
{
  m_reference_count+=1;
  
}

void Directions::Release()
{
  if (m_reference_count<2)
    {
      delete(this);
      
    }
  else
    {
      m_reference_count-=1;  
    }
  

}
