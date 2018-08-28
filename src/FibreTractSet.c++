/************************************************************************

   File: FibreTractSet.c++

   Author: Jennifer Campbell

   Created:

   Revisions:

   Description:

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.

***********************************************************************/
#include <float.h>
#include <vtkByteSwap.h>
#include <unistd.h>

#include "FibreTractSet.h"
#include "otherlibs_include.h"
#include "FibreTract.h"
#include "ImageData.h"
#include "misc.h"


FibreTractSet::FibreTractSet()
{
  //_number_of_tracts=0;
  //_tract_allocation_step=10;
  //_fibre_tract_array=(FibreTract **) malloc(_tract_allocation_step*sizeof(FibreTract *));
  //_fibre_tract_array=(FibreTract **) malloc(_tract_allocation_step*sizeof(FibreTract *));
  //_tracts_allocated=_tract_allocation_step;
  _binary_connectivity_map_init=false;
  //_min_DOTT_map_init=false;
  _cum_min_OTT_map_init=false;
  _tract_scalar_output_mode=no_tract_scalar;
  _point_scalar_output_mode=no_point_scalar;
  _output_geometry_parameters_set=false;
  _save_tracts=true;
  _ASCII=false; // IL added
  _prune_dump=NULL;

}

//used to properly prune branched tracts... not sure about the scalar maps!
void FibreTractSet::ReInitialize()
{
  //DO: reinit the tracts too?  currently freeing them
  int i;

  if (_save_tracts) {
    for (i=0;i<_fibre_tract_array.size();i++) {
      _fibre_tract_array[i]->Release();
    }
  }

  if (_save_tracts && _fibre_tract_array.empty()) { //assume these are all 0 already
    if (_binary_connectivity_map_init) {
      _binary_connectivity_map->SetAllVoxels(0);
    }

    //if (_min_DOTT_map_init)
    //{
    //  _min_DOTT_map->SetAllVoxels(0);
    //}
    if (_cum_min_OTT_map_init) {
      _cum_min_OTT_map->SetAllVoxels(0);
    }
  }
  _fibre_tract_array.clear();
  
  if(_prune_dump) fclose(_prune_dump);
  _prune_dump=NULL;
}

FibreTractSet::FibreTractSet(char *filename)
{
  int i;

  //_tract_allocation_step=10;

  //_fibre_tract_array=(FibreTract **) malloc(_tract_allocation_step*sizeof(FibreTract *));
  //_tracts_allocated=_tract_allocation_step;
  _binary_connectivity_map_init=false;
  //_min_DOTT_map_init=false;
  _cum_min_OTT_map_init=false;
  _tract_scalar_output_mode=no_tract_scalar;
  _point_scalar_output_mode=no_point_scalar;
  _output_geometry_parameters_set=false;
  _prune_dump=NULL;

  //get vtk polydata:
  //get scalars if they exist

  //want: scalars per line, lines, and can get points via indices in the lines.  None of this seems to work

  vtkPolyDataReader* polydata_reader=vtkPolyDataReader::New();
  polydata_reader->SetFileName(filename);
  vtkPolyData *tracts_polydata=polydata_reader->GetOutput();

  //tracts_polydata->BuildLinks(); //no change

  //vtkCellTypes *cell_types;

  //tracts_polydata->GetCellTypes(cell_types);

  //0
  //cout << cell_types->GetNumberOfTypes() << endl;

  //cout << cell_types->GetCellType(0) << endl;


  vtkPoints *tracts_points=tracts_polydata->GetPoints();

  cout << "I can't figure out how to implement this method!  The polydata acts as if empty." << endl;


  //DO: parse the input file without vtk and create structure that way/(i.e. don't use vtk at all: use vtk only for vis.

  //haven't tried: GetPoint(...)

  //float x[3];
  //seg fault here:
  //tracts_points->GetPoint(1,x);

  //cout << "point 1: " << x[0] << x[1] << x[2] << endl;

  //seg fault here:
  //cout << "number of points: " << tracts_points->GetNumberOfPoints() << endl;


  //0:
  //vtkCellArray *tracts_cell_array=tracts_polydata->GetLines();
  //cout << "number of cells: " << tracts_cell_array->GetNumberOfCells() << endl;


  //vtkCellData *tracts_polydata->GetCellData();


  //0:
  //_number_of_tracts=tracts_polydata->GetNumberOfLines();
  //cout << "number of tracts: " <<  _number_of_tracts << endl;


}

//FibreTractSet::FibreTractSet(const FibreTractSet &fibre_tract_set);
//{

//}

FibreTractSet::~FibreTractSet()
{
  int i;

  if (_save_tracts) {

    for (i=0;i<_fibre_tract_array.size();i++) {
      _fibre_tract_array[i]->Release();

    }
  }

  //free(_fibre_tract_array);



  //if (_min_DOTT_map_init)
  //{
  //  _min_DOTT_map->Release();
  //}

  if (_cum_min_OTT_map_init) {
    _cum_min_OTT_map->Release();

  }

  if (_binary_connectivity_map_init) {
    _binary_connectivity_map->Release();

  }

  //DO: plus others: CI

  //Q is something missing here?
  if(_prune_dump)
    fclose(_prune_dump);
}


//note that this method does not automatically update scalar maps
//(do we want to add that?  are there cases where we wouldn't want to?)
void FibreTractSet::AddTract(FibreTract *fibre_tract)
{
  //if (_number_of_tracts==_tracts_allocated) {
  //  this->_ReallocateTractArray();
  //}

  //_fibre_tract_array[_number_of_tracts]=fibre_tract;

  if (_save_tracts) {
    fibre_tract->Retain();
  }
  _fibre_tract_array.push_back(fibre_tract);

  //_number_of_tracts+=1;
}




void FibreTractSet::OutputVTKAscii() //IRL created for ASCII option
{
  _ASCII=true;
}

void FibreTractSet::OutputOBJFile(const char *filename,bool clobber) //TODO: VF finish this
{
  FILE *fp;
  int i,j,k;
  int totpoints=0;
  int offset=0;
  POINT3D point;

  //added by IRL
  INDEX3D index;
  int num_points;
  int j_swap;
  float tmp;

  if (!clobber && !access(filename, F_OK)) {
    printf("Output file %s exists ! Use -clobber!\n",filename);
    return;
  }

  fp=fopen(filename,"w");

  // IL want to write in binary

  for (i=0;i<_fibre_tract_array.size();i++) {
    totpoints+=_fibre_tract_array[i]->GetNumberOfPoints();
  }

  fprintf(fp,"L 0.1 %d\n",totpoints);

  //WRITING cordinates

  for (i=0;i<_fibre_tract_array.size();i++) {
    for (j=0;j<_fibre_tract_array[i]->GetNumberOfPoints();j++) {
      point=_fibre_tract_array[i]->GetPoint(j);
      fprintf(fp,"%f %f %f\n",point.x,point.y,point.z);
    }
  }

  fprintf(fp,"%d\n1 ",_fibre_tract_array.size()); //1 color per line

  for (i=0;i<_fibre_tract_array.size();i++) { //random colors
    fprintf(fp,"%f %f %f 1.0\n",rand()/(double)RAND_MAX,rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
  }

  int cnt=0;

  for (i=0;i<_fibre_tract_array.size();i++) {
    cnt+=_fibre_tract_array[i]->GetNumberOfPoints();
    fprintf(fp,"%d ",cnt);
  }

  fprintf(fp,"\n");

  cnt=0;

  for (i=0;i<_fibre_tract_array.size();i++) {
    //fprintf(fp,"%d\n",cnt);
    for (j=0;j<_fibre_tract_array[i]->GetNumberOfPoints();j++,cnt++)
      fprintf(fp," %d",cnt);

    fprintf(fp,"\n");
  }

  fclose(fp);
}


void FibreTractSet::OutputVTKFile(const char *filename)
{
  FILE *fp;
  int i,j,k;
  int totpoints=0;
  int offset=0;
  POINT3D point;

  //added by IRL
  INDEX3D index;
  int num_points;
  int j_swap;
  float tmp;

  if (!_save_tracts) {
    Error("tracts not saved"); //IRL changed to Error to not get warning "deprecated conversion from string constant to 'char*'"

  }

  fp=fopen(filename,"w");

  if (!_ASCII) {
    // IL want to write in binary
    fprintf(fp,"# vtk DataFile Version 3.0\ntracts\nBINARY\nDATASET POLYDATA\nPOINTS ");

    for (i=0;i<_fibre_tract_array.size();i++) {
      totpoints+=_fibre_tract_array[i]->GetNumberOfPoints();

    }

    fprintf(fp,"%i float\n",totpoints);

    for (i=0;i<_fibre_tract_array.size();i++) {
      //cout << "tract: " << i << " number of points: " << _fibre_tract_array[i]->GetNumberOfPoints() << endl;
      for (j=0;j<_fibre_tract_array[i]->GetNumberOfPoints();j++) {

        point=_fibre_tract_array[i]->GetPoint(j);
        // IL: added index and cout
        //index=_fibre_tract_array[i]->GetIndex(j);
        //cout <<index.x<<" "<<index.y<<" "<<index.z<<endl;
        // IL need to byte swap, no more carriage returns
        vtkByteSwap::Swap4BE(&point.x);
        fwrite(&point.x, sizeof(float),1,fp);
        vtkByteSwap::Swap4BE(&point.y);
        fwrite(&point.y, sizeof(float),1,fp);
        vtkByteSwap::Swap4BE(&point.z);
        fwrite(&point.z, sizeof(float),1,fp);
      }
    }

    fprintf(fp,"\n");

    fprintf(fp,"LINES %i %i\n",_fibre_tract_array.size(),_fibre_tract_array.size()+totpoints);

    for (i=0;i<_fibre_tract_array.size();i++) {
      num_points=_fibre_tract_array[i]->GetNumberOfPoints();
      vtkByteSwap::Swap4BE(&num_points);
      fwrite(&num_points,sizeof(int),1,fp);

      for (j=offset;j<offset+_fibre_tract_array[i]->GetNumberOfPoints();j++) {
        j_swap=j;
        vtkByteSwap::Swap4BE(&j_swap);
        fwrite(&j_swap, sizeof(int),1,fp);
      }

      offset+=_fibre_tract_array[i]->GetNumberOfPoints();

      //fprintf(fp,"\n"); //IL don't need carriage returns

    }
  } else { //IL _ASCII=true (write in ASCII)
    fprintf(fp,"# vtk DataFile Version 3.0\ntracts\nASCII\nDATASET POLYDATA\nPOINTS ");

    for (i=0;i<_fibre_tract_array.size();i++) {
      totpoints+=_fibre_tract_array[i]->GetNumberOfPoints();

    }

    fprintf(fp,"%i float\n",totpoints);

    for (i=0;i<_fibre_tract_array.size();i++) {

      // IL: commented in
      //cout << "tract: " << i << " number of points: " << _fibre_tract_array[i]->GetNumberOfPoints() << endl;
      for (j=0;j<_fibre_tract_array[i]->GetNumberOfPoints();j++) {

        point=_fibre_tract_array[i]->GetPoint(j);
        // IL added index and cout
        //index=_fibre_tract_array[i]->GetIndex(j);
        //cout <<index.x<<" "<<index.y<<" "<<index.z<<endl;
        fprintf(fp,"%f %f %f\n",point.x,point.y,point.z);

      }
    }

    fprintf(fp,"LINES %i %i\n",_fibre_tract_array.size(),_fibre_tract_array.size()+totpoints);

    for (i=0;i<_fibre_tract_array.size();i++) {
      fprintf(fp,"%i ",_fibre_tract_array[i]->GetNumberOfPoints());

      for (j=offset;j<offset+_fibre_tract_array[i]->GetNumberOfPoints();j++) fprintf(fp,"%i ",j);

      offset+=_fibre_tract_array[i]->GetNumberOfPoints();

      fprintf(fp,"\n");

    }

    //only for scalars at each point along tract (not doing this: no reason to at this point.  could be added)
    //fprintf(fp,"POINT_DATA %i\n",totpoints);
    //fprintf(fp,"SCALARS F float 1\n");
    //fprintf(fp,"LOOKUP_TABLE default\n");      
    //for (i=0;i<tracts.numtracts;i++) {
    //  for (j=0;j<tracts.pointsperline[i];j++) {
    //   fprintf(fp,"%f\n",tracts.points[3][j][i]); 
    //  }
    //}
  }

  if (_tract_scalar_output_mode!=no_tract_scalar) { //IL didn't test binary for this option !!!! might need some /n somewhere


    fprintf(fp,"CELL_DATA %i\n",_fibre_tract_array.size());

    switch (_tract_scalar_output_mode) {
    case CI:
      fprintf(fp,"SCALARS connectivity_index float 1\n");
      break;
    case minDOTT:
      fprintf(fp,"SCALARS minDOTT float 1\n");
      break;
    case meanDOTT:
      fprintf(fp,"SCALARS meanDOTT float 1\n");
      break;
    case minOTT:
      fprintf(fp,"SCALARS minOTT float 1\n");
      break;
      //case misc:
      //fprintf(fp,"SCALARS unknown float 1\n");
      //break;
    }

    if (!_ASCII) {
      //step lookup table will be set in FibreTractPlotter
      fprintf(fp,"LOOKUP_TABLE default\n"); //not sure if need this: won't ever actually want the default LUT

      for (i=0;i<_fibre_tract_array.size();i++) {
        switch (_tract_scalar_output_mode) {
          //IL made some mods to write in binary instead of ASCII
        case CI:
          tmp = _fibre_tract_array[i]->GetTractConnectivityIndex();
          vtkByteSwap::Swap4BE(&tmp);
          fwrite(&tmp, sizeof(float),1,fp);
          break;
        case minDOTT:
          tmp = _fibre_tract_array[i]->GetMinDOTT();
          vtkByteSwap::Swap4BE(&tmp);
          fwrite(&tmp, sizeof(float),1,fp);
          break;
        case meanDOTT:
          tmp = _fibre_tract_array[i]->GetMeanDOTT();
          vtkByteSwap::Swap4BE(&tmp);
          fwrite(&tmp, sizeof(float),1,fp);
          break;
        case minOTT:
          tmp = _fibre_tract_array[i]->GetMinOTT();
          vtkByteSwap::Swap4BE(&tmp);
          fwrite(&tmp, sizeof(float),1,fp);
          break;
        }
      }
    } else { //IL write in ASCII
      //step lookup table will be set in FibreTractPlotter
      fprintf(fp,"LOOKUP_TABLE default\n"); //not sure if need this: won't ever actually want the default LUT

      for (i=0;i<_fibre_tract_array.size();i++) {
        switch (_tract_scalar_output_mode) {
        case CI:
          fprintf(fp,"%f\n",_fibre_tract_array[i]->GetTractConnectivityIndex());
          break;
        case minDOTT:
          fprintf(fp,"%f\n",_fibre_tract_array[i]->GetMinDOTT());
          break;
        case meanDOTT:
          fprintf(fp,"%f\n",_fibre_tract_array[i]->GetMeanDOTT());
          break;
        case minOTT:
          fprintf(fp,"%f\n",_fibre_tract_array[i]->GetMinOTT());
          break;
        }
      }
    }

  }

  if (_point_scalar_output_mode!=no_point_scalar) {


    fprintf(fp,"POINT_DATA %i\n",totpoints);

    switch (_point_scalar_output_mode) {
    case cumulative_CI:
      fprintf(fp,"SCALARS cumulative_connectivity_index float 1\n");
      break;
    case DOTT:
      fprintf(fp,"SCALARS DOTT float 1\n");
      break;
    case cumulative_min_DOTT:
      fprintf(fp,"SCALARS cumulative_min_DOTT float 1\n");
      break;
    case OTT:
      fprintf(fp,"SCALARS OTT float 1\n");
      break;
    case cumulative_min_OTT:
      fprintf(fp,"SCALARS cumulative_min_OTT float 1\n");
      break;
      /*
      case misc_point_scalar:
      fprintf(fp,"SCALARS unknown float 1\n");
      break;
      */
    }

    if (!_ASCII) {
      //step lookup table will be set in FibreTractPlotter
      fprintf(fp,"LOOKUP_TABLE default\n"); //not sure if need this: won't ever actually want the default LUT

      for (i=0;i<_fibre_tract_array.size();i++) {
        for (j=0;j<_fibre_tract_array[i]->GetNumberOfPoints();j++) {
          //IL made some mods to write in binary instead of ASCII
          switch (_point_scalar_output_mode) {
          case cumulative_CI:
            tmp =_fibre_tract_array[i]->GetCumulativeCI(j);
            vtkByteSwap::Swap4BE(&tmp);
            fwrite(&tmp, sizeof(float),1,fp);
            break;
          case DOTT:
            tmp = _fibre_tract_array[i]->GetDOTT(j);
            vtkByteSwap::Swap4BE(&tmp);
            fwrite(&tmp, sizeof(float),1,fp);
            break;

          case cumulative_min_DOTT:
            tmp = _fibre_tract_array[i]->GetCumulativeMinDOTT(j);
            vtkByteSwap::Swap4BE(&tmp);
            fwrite(&tmp, sizeof(float),1,fp);
            break;
          case OTT:
            tmp = _fibre_tract_array[i]->GetDOTT(j);
            vtkByteSwap::Swap4BE(&tmp);
            fwrite(&tmp, sizeof(float),1,fp);
            break;

          case cumulative_min_OTT:
            tmp = _fibre_tract_array[i]->GetCumulativeMinOTT(j);
            vtkByteSwap::Swap4BE(&tmp);
            fwrite(&tmp, sizeof(float),1,fp);
            break;
          }
        }
      }
    } else { //IL write in ASCII
      //step lookup table will be set in FibreTractPlotter
      fprintf(fp,"LOOKUP_TABLE default\n"); //not sure if need this: won't ever actually want the default LUT

      for (i=0;i<_fibre_tract_array.size();i++) {
        for (j=0;j<_fibre_tract_array[i]->GetNumberOfPoints();j++) {

          switch (_point_scalar_output_mode) {
          case cumulative_CI:
            fprintf(fp,"%f\n",_fibre_tract_array[i]->GetCumulativeCI(j));
            break;
          case DOTT:
            fprintf(fp,"%f\n",_fibre_tract_array[i]->GetDOTT(j));
            break;

          case cumulative_min_DOTT:
            fprintf(fp,"%f\n",_fibre_tract_array[i]->GetCumulativeMinDOTT(j));
            break;
          case OTT:
            fprintf(fp,"%f\n",_fibre_tract_array[i]->GetOTT(j));
            break;

          case cumulative_min_OTT:
            fprintf(fp,"%f\n",_fibre_tract_array[i]->GetCumulativeMinOTT(j));
            break;
          }
        }
      }
    }

  }



  fclose(fp);


}


void FibreTractSet::OutputVTKFile(const char *filename, bool clobber)
{
  bool may_output=true;
  char use_new_filename;
  char new_filename[1000];

  //makes the buffer bigger:
  strcpy(new_filename,filename);

  while (may_output) {

    if(clobber || access(new_filename,F_OK))
    {
      this->OutputVTKFile(new_filename);
      break;
    }   else { //prompt user for new filename
        cout << "OutputVTKFile:"<<new_filename<<" exists!"<<std::endl;
      
        cout << "Do you wish to use an alternate file name? [y] [n]"<<std::endl;
        cin  >> use_new_filename;

      if (use_new_filename=='y') {
        cout << "new file name:"<<std::endl;
        cin >> new_filename;
      } else {
        cout<<"Overwriting..."<<std::endl; //IRL I WANT to save!
        clobber=true;
      }
    }
  }

}



void FibreTractSet::SetImageDataTemplate(ImageData *image_data)
{
  _output_geometry_parameters.x_size=image_data->GetXSize();
  _output_geometry_parameters.y_size=image_data->GetYSize();
  _output_geometry_parameters.z_size=image_data->GetZSize();
  _output_geometry_parameters.x_start=image_data->GetXStart();
  _output_geometry_parameters.y_start=image_data->GetYStart();
  _output_geometry_parameters.z_start=image_data->GetZStart();
  _output_geometry_parameters.x_step=image_data->GetXStep();
  _output_geometry_parameters.y_step=image_data->GetYStep();
  _output_geometry_parameters.z_step=image_data->GetZStep();
  _output_geometry_parameters.x_direction_cosine=image_data->GetXDirectionCosine();
  _output_geometry_parameters.y_direction_cosine=image_data->GetYDirectionCosine();
  _output_geometry_parameters.z_direction_cosine=image_data->GetZDirectionCosine();

  _output_geometry_parameters_set=true;

}

void FibreTractSet::SetImageDataTemplate(const char *image_filename)
{
  ImageData *image_data = new ImageData(image_filename);

  SetImageDataTemplate(image_data);

  image_data->Release();
}

void FibreTractSet::SetImageDataParameters(int x_size, int y_size, int z_size,
    float x_start, float y_start, float z_start,
    float x_step, float y_step, float z_step,
    const VECTOR3D& x_direction_cosine,
    const VECTOR3D& y_direction_cosine,
    const VECTOR3D& z_direction_cosine)
{
  _output_geometry_parameters.x_size=x_size;
  _output_geometry_parameters.y_size=y_size;
  _output_geometry_parameters.z_size=z_size;
  _output_geometry_parameters.x_start=x_start;
  _output_geometry_parameters.y_start=y_start;
  _output_geometry_parameters.z_start=z_start;

  _output_geometry_parameters.x_step=x_step;
  _output_geometry_parameters.y_step=y_step;
  _output_geometry_parameters.z_step=z_step;

  _output_geometry_parameters.x_direction_cosine=x_direction_cosine;
  _output_geometry_parameters.y_direction_cosine=y_direction_cosine;
  _output_geometry_parameters.z_direction_cosine=z_direction_cosine;

  _output_geometry_parameters_set=true;
}

void FibreTractSet::SetImageDataParameters(IMAGE_GEOM image_geom)
{
  _output_geometry_parameters.x_size=image_geom.x_size;
  _output_geometry_parameters.y_size=image_geom.y_size;
  _output_geometry_parameters.z_size=image_geom.z_size;
  _output_geometry_parameters.x_start=image_geom.x_start;
  _output_geometry_parameters.y_start=image_geom.y_start;
  _output_geometry_parameters.z_start=image_geom.z_start;
  _output_geometry_parameters.x_step=image_geom.x_step;
  _output_geometry_parameters.y_step=image_geom.y_step;
  _output_geometry_parameters.z_step=image_geom.z_step;
  _output_geometry_parameters.x_direction_cosine=image_geom.x_direction_cosine;
  _output_geometry_parameters.y_direction_cosine=image_geom.y_direction_cosine;
  _output_geometry_parameters.z_direction_cosine=image_geom.z_direction_cosine;
  _output_geometry_parameters_set=true;
}


void FibreTractSet::InitBinaryConnectivityMap()
{

  if (!_output_geometry_parameters_set) {
    Error("output geometry parameters not set");
  }


  if (!_binary_connectivity_map_init) {

    _binary_connectivity_map=new
    ImageData(Image3D_vanilla,
              _output_geometry_parameters.x_size,_output_geometry_parameters.y_size,_output_geometry_parameters.z_size,
              _output_geometry_parameters.x_start,_output_geometry_parameters.y_start,_output_geometry_parameters.z_start,
              _output_geometry_parameters.x_step,_output_geometry_parameters.y_step,_output_geometry_parameters.z_step,
              _output_geometry_parameters.x_direction_cosine,_output_geometry_parameters.y_direction_cosine,_output_geometry_parameters.z_direction_cosine);

    _binary_connectivity_map->SetAllVoxels(0);


    _binary_connectivity_map_init=true;

  }
}

//this can be run without init to 0. (faster)
void FibreTractSet::ComputeBinaryConnectivityMap()
{
  float binary_conn;
  int x,y,z;
  int i;

  if (!_binary_connectivity_map_init) {
    if (!_output_geometry_parameters_set) {
      Error("output geometry parameters not set");
    }

    _binary_connectivity_map=new ImageData(Image3D_vanilla,_output_geometry_parameters.x_size,_output_geometry_parameters.y_size,_output_geometry_parameters.z_size,_output_geometry_parameters.x_start,_output_geometry_parameters.y_start,_output_geometry_parameters.z_start,_output_geometry_parameters.x_step,_output_geometry_parameters.y_step,_output_geometry_parameters.z_step,_output_geometry_parameters.x_direction_cosine,_output_geometry_parameters.y_direction_cosine,_output_geometry_parameters.z_direction_cosine);

    _binary_connectivity_map_init=true;

  }

  //faster: go through all tracts and round points (problems with points on the edge though...

  for (x=0;x<_binary_connectivity_map->GetXSize();x++) {
    for (y=0;y<_binary_connectivity_map->GetYSize();y++) {
      for (z=0;z<_binary_connectivity_map->GetZSize();z++) {
        binary_conn=0;

        for (i=0;i<_fibre_tract_array.size() && binary_conn<FLT_EPSILON;i++) {
          if (this->GetTract(i)->GoesThroughVoxel(_binary_connectivity_map,x,y,z)) {
            binary_conn=1;


          }

        }

        _binary_connectivity_map->SetValue(x,y,z,binary_conn);

      }
    }
  }



}



ImageData *FibreTractSet::GetBinaryConnectivityMap()
{
  return _binary_connectivity_map;



}




void FibreTractSet::InitCumMinOTTMap()
{
  if (!_output_geometry_parameters_set) {
    Error("output geometry parameters not set");
  }

  if (!_cum_min_OTT_map_init) {
    _cum_min_OTT_map=new ImageData(Image3D_vanilla,_output_geometry_parameters.x_size,_output_geometry_parameters.y_size,_output_geometry_parameters.z_size,_output_geometry_parameters.x_start,_output_geometry_parameters.y_start,_output_geometry_parameters.z_start,_output_geometry_parameters.x_step,_output_geometry_parameters.y_step,_output_geometry_parameters.z_step,_output_geometry_parameters.x_direction_cosine,_output_geometry_parameters.y_direction_cosine,_output_geometry_parameters.z_direction_cosine);

    _cum_min_OTT_map->SetAllVoxels(0);


    _cum_min_OTT_map_init=true;

  }
}

/*
//not sure why this is commented out?
void FibreTractSet::ComputeCumMinOTTMap()
{

  int i,j;
  INDEX3D index;
  FibreTract *tract;

  if (!_cum_min_OTT_map_init)
    {
      this->InitCumMinOTTMap()
    }


  //for each tract,
  for (i=0;i<_number_of_tracts)
    {
      tract=this->GetTract(i);

      //tract->AssignCumulativeMinOTTs();
      for (j=0;j<tract->GetNumberOfPoints();j++)
 {
   //step through tract and take the CMOTTs from it and assign to the map at that index, if greater than current value.
   index=tract->GetPointInfo(j).index;
   if (tract->GetPointInfo(j).cumulative_min_OTT>_cum_min_OTT_map->GetValue(index.x,index.y,index.z))
     {
       _cum_min_OTT_map->SetValue(index.x,index.y,index.z,tract->GetPointInfo(j).cumulative_min_OTT);
     }
 }
    }


}
*/


ImageData *FibreTractSet::GetCumMinOTTMap()
{
  //maybe check if init...
  return _cum_min_OTT_map;

}

FibreTractSet *FibreTractSet::Prune(float threshold)
{
  //keep tracts with scalar index above and including threshold.
  Error("not implemented");
}


//VF: this function takes 90% of execution time, it should be optimized
FibreTractSet *FibreTractSet::PruneSameIndices(bool compare_last_only)
{
  //get rid of any tracts that have exactly the same list of voxels passed through; keep the most "likely"
  
  if (!_save_tracts) {
    Error("Can't prune: tracts are not saved");//IRL changed to Error to get rid of warning

  }
  /*
  if(!_prune_dump)
  {
    char tmp[128];
    sprintf(tmp,"prune_dump_%d.csv",(unsigned long)getpid());
    _prune_dump=fopen(tmp,"a");
  }*/
  int i,j,k;

  int lo,sh,t,q,start_index; //IL: indices to loop over short and long tracts

  FibreTractSet *pruned_fibre_tract_set = new FibreTractSet();
  FibreTract *short_tract;
  FibreTract *long_tract;
  FibreTract *new_tract;

  bool stop=false;
  bool add_last_tract_at_end=true;
  bool found=false; //IL: becomes true if we find a subset
  bool current_is_long=false; //IL: becomes true when current tract has already replaced a subset
  float x,y,z;

  if (compare_last_only) {
    
    for (i=0;i<(GetNumberOfTracts()-1) && !stop;i++) 
    {

      // IL: get length of tracts to compare
      if (GetTract(GetNumberOfTracts()-1)->GetNumberOfPoints() <= GetTract(i)->GetNumberOfPoints()) 
      {
        short_tract=GetTract(this->GetNumberOfTracts()-1);
        long_tract=GetTract(i);
      }  else {
        short_tract=GetTract(i);
        long_tract=GetTract(GetNumberOfTracts()-1);
      }

      sh = short_tract->GetNumberOfPoints();
      lo = long_tract->GetNumberOfPoints();
      
      //fprintf(_prune_dump,"%d,%d\n",lo,sh);

      found=false;

      for (t=0; t<(lo-sh) && !found; t++) { // loop over long tract , check if short is subset
        
        for (j=0; (t+j)<lo  && j<sh && short_tract->GetIndex(j)==long_tract->GetIndex(t+j);j++) 
        {

          //Ilana debug begin
          /*cout<<"long tract length: "<<lo<<"short tract length:"<<sh<<endl;
          cout <<"j t "<<j<<" "<<t<<endl;
          cout <<"x index: short  long"<<short_tract->GetIndex(j).x<<" "<<long_tract->GetIndex(t+j).x<<endl;
          cout <<"x index: short  long"<<short_tract->GetIndex(j).y<<" "<<long_tract->GetIndex(t+j).y<<endl;
          cout <<"x index: short  long"<<short_tract->GetIndex(j).z<<" "<<long_tract->GetIndex(t+j).z<<"\n\n"<<endl;*/
        }

        if (j>=sh) { // short tract is a subset of long
          add_last_tract_at_end=false;
          found=true;
          start_index = t; // store index at which short tract starts in long tract
        }
      }  // t loop for long_tract

      if (found) {
        if (GetTract(GetNumberOfTracts()-1)->GetNumberOfPoints() <= GetTract(i)->GetNumberOfPoints()) {
          //short_tract=this->GetTract(this->GetNumberOfTracts()-1);
          //long_tract=this->GetTract(i);
          stop=true;
        } else {
          //short_tract=this->GetTract(i);
          //long_tract=this->GetTract(this->GetNumberOfTracts()-1);
        }

        if (!current_is_long) { //IL: have to inbed this if statement or else segmentation fault
          pruned_fibre_tract_set->AddTract(long_tract); // only add this tract if it was not already added

          //cout<<"added long tract "<<_number_of_tracts <<" with "<<long_tract->GetNumberOfPoints()<<" at place of tract "<<i<<" with "<<short_tract->GetNumberOfPoints()<<endl;
          //for (k=0;k<long_tract->GetNumberOfPoints();k++)
          //{
          //    cout<<long_tract->GetIndex(k).x<<" "<<long_tract->GetIndex(k).y<<" "<<long_tract->GetIndex(k).z<<endl;
          //}

          if (GetTract(GetNumberOfTracts()-1)->GetNumberOfPoints()>GetTract(i)->GetNumberOfPoints()) {
            current_is_long=true; // we only want to add new tract once
          }
        }
      } else  {
        //else, they are not the same://IL: keep original tract in set (i.e. just transfer to output set)
        pruned_fibre_tract_set->AddTract(GetTract(i));
        //cout<<"tract not found added tract "<<i<<" with "<<this->GetTract(i)->GetNumberOfPoints()<<endl;
      }
    } // loop i over number of tracts

    if(stop) {
      //add the rest of the old tracts:
      for (i=i;i<(_fibre_tract_array.size()-1);i++) {
        pruned_fibre_tract_set->AddTract(GetTract(i));
        //cout<<"just added old tract "<<i<<" with "<<this->GetTract(i)->GetNumberOfPoints()<<endl;
      }
    }

    if (add_last_tract_at_end) {
      //add the (new) last one:
      pruned_fibre_tract_set->AddTract(GetTract(GetNumberOfTracts()-1));
      //cout <<"just added tract at end "<<this->GetNumberOfTracts()-1<<" with "<<this->GetTract(this->GetNumberOfTracts()-1)->GetNumberOfPoints()<<endl;
    }


    //cout <<   pruned_fibre_tract_set->GetNumberOfTracts() << endl;
  } else {
    //compare all tracts to all tracts: not done
    Error("not implemented: must prune comparing last only.");

  }



  //all of this is basically a deep copy with everything except the tracts:

  if (_output_geometry_parameters_set) {
    pruned_fibre_tract_set->SetImageDataParameters(_output_geometry_parameters);

  }

  //assign the image datas: note that binary conn map is (in FibreTracking) being assigned on the fly so Prune should not change it


  //if (_min_DOTT_map_init)
  //{
  //  pruned_fibre_tract_set->SetMinDOTTMap(_min_DOTT_map);
  //}

  if (_cum_min_OTT_map_init) {
    pruned_fibre_tract_set->SetCumMinOTTMap(_cum_min_OTT_map);
  }


  if (_binary_connectivity_map_init) {
    pruned_fibre_tract_set->SetBinaryConnectivityMap(_binary_connectivity_map);
  }


  pruned_fibre_tract_set->SetTractScalarOutputMode(_tract_scalar_output_mode);

  pruned_fibre_tract_set->SetPointScalarOutputMode(_point_scalar_output_mode);


  //end deep copy

  return pruned_fibre_tract_set;
}


//this uses the tract indices to select only the tracts that go through the roi:
FibreTractSet *FibreTractSet::Prune(ImageData *roi, bool exclude, bool compare_last_only)
{
  int i,k;
  int x,y,z;
  INDEX3D index;
  bool tract_added;


  if (!_save_tracts) {
    Error("Can't prune: tracts not saved");

  }


  FibreTractSet *pruned_fibre_tract_set = new FibreTractSet();

  if (_output_geometry_parameters_set) {
    pruned_fibre_tract_set->SetImageDataParameters(_output_geometry_parameters);

  }

  //if (_min_DOTT_map_init)
  //{
  //  pruned_fibre_tract_set->InitMinDOTTMap();
  //}


  if (_cum_min_OTT_map_init) {
    pruned_fibre_tract_set->InitCumMinOTTMap();
  }


  if (_binary_connectivity_map_init) {
    pruned_fibre_tract_set->InitBinaryConnectivityMap();
  }


  if (!compare_last_only) {
    if (exclude) {
      for (i=0;i<_fibre_tract_array.size();i++) {
        tract_added=true;

        for (k=0;k<this->GetTract(i)->GetNumberOfPoints() && tract_added ;k++) {

          if (roi->GetValue(this->GetTract(i)->GetIndex(k).x,this->GetTract(i)->GetIndex(k).y,this->GetTract(i)->GetIndex(k).z)>FLT_EPSILON) {


            tract_added=false;


          }
        }

        if (tract_added) {
          pruned_fibre_tract_set->AddTract(this->GetTract(i));
        }


      }
    } else {
      for (i=0;i<_fibre_tract_array.size();i++) {
        tract_added=false;

        for (k=0;k<this->GetTract(i)->GetNumberOfPoints() && !tract_added ;k++) {

          if (roi->GetValue(this->GetTract(i)->GetIndex(k).x,this->GetTract(i)->GetIndex(k).y,this->GetTract(i)->GetIndex(k).z)>FLT_EPSILON) {

            pruned_fibre_tract_set->AddTract(this->GetTract(i));
            tract_added=true;


          }
        }

      }
    }
  } else { //compare_last_only
    if (exclude) {
      for (i=_fibre_tract_array.size()-1;i<_fibre_tract_array.size();i++) {
        tract_added=true;

        for (k=0;k<this->GetTract(i)->GetNumberOfPoints() && tract_added ;k++) {

          if (roi->GetValue(this->GetTract(i)->GetIndex(k).x,this->GetTract(i)->GetIndex(k).y,this->GetTract(i)->GetIndex(k).z)>FLT_EPSILON) {


            tract_added=false;


          }
        }

        if (tract_added) {
          pruned_fibre_tract_set->AddTract(this->GetTract(i));
        }


      }
    } else {
      for (i=_fibre_tract_array.size()-1;i<_fibre_tract_array.size();i++) {
        tract_added=false;

        for (k=0;k<this->GetTract(i)->GetNumberOfPoints() && !tract_added ;k++) {

          if (roi->GetValue(this->GetTract(i)->GetIndex(k).x,this->GetTract(i)->GetIndex(k).y,this->GetTract(i)->GetIndex(k).z)>FLT_EPSILON) {

            pruned_fibre_tract_set->AddTract(this->GetTract(i));
            tract_added=true;


          }
        }

      }
    }



  } //compare_last_only









  if (_output_geometry_parameters_set && (_binary_connectivity_map_init || _cum_min_OTT_map_init)) {

    for (i=0;i<pruned_fibre_tract_set->GetNumberOfTracts();i++) {
      for (k=0;k<pruned_fibre_tract_set->GetTract(i)->GetNumberOfPoints();k++) {
        index=pruned_fibre_tract_set->GetTract(i)->GetIndex(k);

        if (_binary_connectivity_map_init) {
          pruned_fibre_tract_set->GetBinaryConnectivityMap()->SetValue(index.x,index.y,index.z,this->GetBinaryConnectivityMap()->GetValue(index.x,index.y,index.z));
        }

        //if (_min_DOTT_map_init)
        //{
        //  pruned_fibre_tract_set->GetMinDOTTMap()->SetValue(index.x,index.y,index.z,this->GetMinDOTTMap()->GetValue(index.x,index.y,index.z));
        //}
        if (_cum_min_OTT_map_init) {
          pruned_fibre_tract_set->GetCumMinOTTMap()->SetValue(index.x,index.y,index.z,this->GetCumMinOTTMap()->GetValue(index.x,index.y,index.z));
        }


      }
    }
  }




  pruned_fibre_tract_set->SetTractScalarOutputMode(_tract_scalar_output_mode);

  pruned_fibre_tract_set->SetPointScalarOutputMode(_point_scalar_output_mode);



  //end deep copy

  return pruned_fibre_tract_set;


}


void FibreTractSet::SetTractScalarOutputMode(TRACT_SCALAR_OUTPUT_MODE tract_scalar_output_mode)
{
  _tract_scalar_output_mode=tract_scalar_output_mode;

}

void FibreTractSet::SetPointScalarOutputMode(POINT_SCALAR_OUTPUT_MODE point_scalar_output_mode)
{
  _point_scalar_output_mode=point_scalar_output_mode;

}


/*

void FibreTractSet::_ReallocateTractArray()
{

  _fibre_tract_array=(FibreTract **) realloc(_fibre_tract_array,(_tracts_allocated+_tract_allocation_step)*sizeof(FibreTract *));
  _tracts_allocated+=_tract_allocation_step;

  _tract_allocation_step*=2;

}*/

void FibreTractSet::SetCumMinOTTMap(ImageData *cum_min_OTT_map)
{
  _cum_min_OTT_map_init=true; //assuming the argument is allocated!
  _cum_min_OTT_map=cum_min_OTT_map;
  _cum_min_OTT_map->Retain();
}


void FibreTractSet::SetBinaryConnectivityMap(ImageData *binary_connectivity_map)
{
  _binary_connectivity_map_init=true;
  _binary_connectivity_map=binary_connectivity_map;
  _binary_connectivity_map->Retain();
}

void FibreTractSet::TractsOff()
{
  _save_tracts=false;

}

bool FibreTractSet::SaveTractsOn()
{
  return _save_tracts;


}



//this removes a tract but assumes scalars stay as they were: it is just for removing a tract that is a subset of other tracts in this set.
void FibreTractSet::RemoveTract(int n)
{
  int i;

  if (_save_tracts) 
  {
    _fibre_tract_array[n]->Release();

    /*
    for (i=n;i<_fibre_tract_array.size()-1;i++) {
      _fibre_tract_array[n]=_fibre_tract_array[n+1];
    }
    _fibre_tract_array.pop_back();*/
    _fibre_tract_array.erase(_fibre_tract_array.begin()+n,_fibre_tract_array.begin()+n+1);
  }
}






























