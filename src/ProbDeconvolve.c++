	/************************************************************************
	
	   File: ProbDeconvolve.c++
	
	   Author: Parya Momayyez & Jennifer Campbell
	
	   Created: Aug 4th, 2008
	
	
	   Description:
	
	   Copyright (c) 2008, 2011 by Parya Momayyez & Jennifer Campbell,
	   McConnell Brain Imaging Center and Centre for Intelligent Machines,
	   McGill University, Montreal, QC.
	
	 ***********************************************************************/
	
	#ifndef __jc_tools_h
	#include "jc_tools.h"
	#endif
	
	#ifndef __otherlibs_include_h
	#include "otherlibs_include.h"
	#endif
	
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_sf.h>
	#include <gsl/gsl_statistics.h>
	#include <gsl/gsl_blas.h>
	
	/*
	

!This program has memory trouble if the mask is sufficiently large; may need to do in chunks in that case

this program and all of the HAR stuff assume data acquired on one hemisphere.
	
	*/
	

	
//NOTE!  this doesn't calculate maxima and sigma anymore; the bootstrap samples instead
void CalculateMaximaAndSigma(char *input_name, char *mask_name, int order, char *output_name, int iterations, float lambda1, float lambda2, float So, bool clobber, char *history)
	{
	  // Declaration of variables for the residual bootstrapping and deconvolution part
	
	  float maskVal=0.5; //FLT_EPSILON; //parya had input this, which makes more sense: can add to mask input section 
		int i,j,k,x,y,z,n,m;
		int seed;
	
		int t;
	
		seed = (int) time (NULL);
		srand48(seed);
	
		DeconvolutionAlgorithm algorithm = FORECAST;
		Deconvolver **Deconvolver_array = (Deconvolver **) malloc(iterations*sizeof(Deconvolver *));
	
		ImageData *DWI_image_data = new ImageData(input_name);
		ImageData *mask_image_data = new ImageData(mask_name);
		ImageData *residuals;
		ImageData **DWI_data_array = (ImageData **) malloc(iterations*sizeof(ImageData *));
		ImageData **Deconv_data_array = (ImageData **) malloc(iterations*sizeof(ImageData *));

		ImageData *tmp_image;
		//ImageData *response = new ImageData(input_name);
	
		//ImageData *response_voxel = new ImageData(response_mask_name);
	
	
		//ImageData *lambda1;
		//ImageData *lambda2;
	
		DWI_data_array[0] = DWI_image_data->CreateDWIODF();
		
	
		if (DWI_data_array[0]->GetNSize()<50)
		  {
		    tmp_image=Create100DirInput(DWI_data_array[0]);
		    DWI_data_array[0]->Release();
		    DWI_data_array[0]=tmp_image;
		  }
	
		int xsize,ysize,zsize;
		xsize=DWI_image_data->GetXSize();
		ysize=DWI_image_data->GetYSize();
		zsize=DWI_image_data->GetZSize();
	
		SphericalHarmonics *sphere_harm = new SphericalHarmonics();
		sphere_harm->SetOrder(order);
		ImageData *SH_fit = new ImageData(ODFdata,DWI_data_array[0]->GetImageGeom());
		residuals = new ImageData(ODFdata,DWI_data_array[0]->GetImageGeom());
	
	
		float *lev;
			
		int res_index1,res_index2;
		float temp_res;

		//full_avg_odf for fanning:

		//hardcoded not to do the actual distributions; avg the odf:
		int justpdf=false;

		ImageData *full_avg_odf;

		/*
		if (algorithm == FORECAST)//==only option right now
		{
			lambda1 = new ImageData(lambda1_name);
			lambda2 = new ImageData(lambda2_name);
		}
		else //Tournier algorithm: not implemented
		{
	
	
	
		}
		*/

		float max_for_norm;

		//float value_for_vol_one=pow((DWI_data_array[0]->GetDirections()->GetNumberOfDirections()*6.0)/(4.0*M_PI),1.0/3.0);

		//for vis for petrides 01, we need 1.0 to be consistent with the vectors
		//uncomment for new datasets:
		//float value_for_vol_one=pow((DWI_data_array[0]->GetDirections()->GetNumberOfDirections()*6.0)/(4.0*M_PI),1.0/3.0);
		//comment for new datasets:
		float value_for_vol_one=1.0;


		for (x=0; x<xsize; x++)
		{
			for (y=0; y<ysize; y++)
			{
				for (z=0; z<zsize; z++)
				{
					if (mask_image_data->GetValue(x,y,z) > maskVal)
					{

					  
						sphere_harm->CalculateSHFit(DWI_data_array[0], SH_fit, x, y, z);

						
						
						lev=(float *) malloc(DWI_data_array[0]->GetNSize()*sizeof(float));
						sphere_harm->GetLeverageTerms(lev,x,y,z,DWI_data_array[0]->GetNSize());
					


						/*test
						for (n=0; n<DWI_data_array[0]->GetNSize(); n++)
						  {
						    lev[n]=0.0;
						  }
						*/
						for (n=0; n<DWI_data_array[0]->GetNSize(); n++)
						{
						  
						 
						  
						  residuals->SetValue(x, y, z, n, (SH_fit->GetValue(x, y, z, n) - DWI_data_array[0]->GetValue(x, y, z, n))/sqrt(1.0-lev[n]));
       
						 
						  
						  //cout << lev[n] << " " << residuals->GetValue(x, y, z, n)<< " " << DWI_data_array[0]->GetValue(x, y, z, n) << " " << SH_fit->GetValue(x, y, z, n) << endl;
	
					
	
						}
						free(lev);
					}
				}
			}
		}
	
		
	
		delete(sphere_harm);
	
		//from parya's code:
	
		// Declaration of variables for the probabilistic part
		
		//for now there are no labels
		bool labels_exist = false;

		float value, dotProd, labelCount, maxDot[2], dot_matrix[3][3];
		int xyznumdirs, ind, sInd[3], mInd[3];
		int labelNum, orientNum;

		if (labels_exist)
		{
			labelNum = 4;
			orientNum = 7;
		}
		else
		{
			labelNum = 3;
			orientNum = 6;
		}

		bool labelsCheck[orientNum];
		int num_of_dir[labelNum];
	
		num_of_dir[0] = 1;
		num_of_dir[1] = 2;
		num_of_dir[2] = 3;
		if (labels_exist)
		{
			num_of_dir[3] = 1;
		}
	
		VECTOR3D dir, dir2, mdir[3];
		Directions *Dirs;
	
		// An array of six/seven imageData objects to save the histogram corresponding to each of the maximum for all possible configurations.
		// Single (1 maximum) + Double Crossing (2 maxima) + Triple Crossing (3 maxima) (+fanning) = 6(+1) maxima in total
		// Description of fourth dimension: 6(Label occurrence rate + Maxima occurrence rate + 95% confidence interval + Mean direction) + number of axonal bins
	
		ImageData **histogram_image_array = (ImageData **) malloc(orientNum * sizeof(ImageData *));
		IMAGE_GEOM img_geom = DWI_data_array[0]->GetImageGeom();
		img_geom.n_size = DWI_data_array[0]->GetDirections()->GetNumberOfDirections()+6;





		for(i=0; i<orientNum; i++)
		{
			histogram_image_array[i] = new ImageData(ODFdata, img_geom);
			histogram_image_array[i]->SetDirections(DWI_data_array[0]->GetDirections());
	
			histogram_image_array[i]->SetAllVoxels(0);
		}
	

		//the full ODF image:

		full_avg_odf = new ImageData(ODFdata, DWI_data_array[0]->GetImageGeom());
		full_avg_odf->SetDirections(DWI_data_array[0]->GetDirections());
		full_avg_odf->SetAllVoxels(0);

		float *dataArray = (float *) malloc(orientNum * xsize * ysize * zsize * img_geom.n_size * sizeof(float));
		float volume[6] = {0, 0, 0, 0, 0, 0};
		int dataInd, offInd;
	
		for (i=0; i<iterations; i++)
		{
		  //cout << "\r" << i; //JC
			if (i>0)
			{
				DWI_data_array[i]=new ImageData(DWIdata, DWI_data_array[0]->GetImageGeom());
				DWI_data_array[i]->SetDiffusionParameters(DWI_data_array[0]->GetDiffusionParameters());
			}
	
			Deconvolver_array[i]=new Deconvolver();
			Deconvolver_array[i]->SetOrder(order);
			Deconvolver_array[i]->SetDeconvolutionAlgorithm(algorithm);
	
			
			if (algorithm==FORECAST)//==only option
			{
				//Deconvolver_array[i]->SetResponseFunction(response);
				//response->Retain();
			  //cout << lambda1 << " " << lambda2 << " " << So << endl;
			  Deconvolver_array[i]->SetLambdas(lambda1,lambda2,So);
			}
			else //Tournier algorithm
			{
	
			}
			
			if (i>0)
			{
			    for (x=0; x<xsize; x++)
			    {
					for (y=0; y<ysize; y++)
					{
						for (z=0; z<zsize; z++)
						{
							if (mask_image_data->GetValue(x,y,z) > 0.5)
							{
								for (n=0; n<DWI_data_array[0]->GetNSize(); n++)
								{
									res_index1 = (int) ((float) DWI_data_array[0]->GetNSize()) * ((float) ((double) rand() / ( (double)(RAND_MAX) + (double)(1.0)) ));
									res_index1 = int_min(res_index1, DWI_data_array[0]->GetNSize()-1);
									//cout << SH_fit->GetValue(x, y, z, n) << endl;
									if (SH_fit->GetValue(x, y, z, n) + residuals->GetValue(x, y, z, res_index1)<0)
									{
										//cout << "! negative value in bootstrap data: SH_fit " <<  SH_fit->GetValue(x, y, z, n) << " residual " << residuals->GetValue(x, y, z, res_index1) << endl;
										DWI_data_array[i]->SetValue(x, y, z, n, 0.0);//not recommended to use this processing if this occurs frequently
									}
									else
									{
										DWI_data_array[i]->SetValue(x, y, z, n, SH_fit->GetValue(x, y, z, n) + residuals->GetValue(x, y, z, res_index1));
										//cout << SH_fit->GetValue(x, y, z, n)  << " " << residuals->GetValue(x, y, z, res_index1) << endl;
									}
								}
							}
						}
					}
				}
			}
			Deconvolver_array[i]->SetInput(DWI_data_array[i]);
	
			for (x=0; x<xsize; x++)
			{
				for (y=0; y<ysize; y++)
				{
					for (z=0; z<zsize; z++)
					{
						if (mask_image_data->GetValue(x,y,z) > 0.5)
						{
						  //cout << "deconvolve \n";
							
							Deconvolver_array[i]->Deconvolve(x,y,z);
						}
					}
				}
			}

			Deconv_data_array[i] = Deconvolver_array[i]->GetOutput();

		

		



			delete(Deconvolver_array[i]);
			if (i>0)
			{
				DWI_data_array[i]->Release();
			}
	
	
			for (x=0; x<xsize; x++)
			{
				for (y=0; y<ysize; y++)
				{
					for (z=0; z<zsize; z++)
					{
						if (mask_image_data->GetValue(x,y,z) > maskVal)
						{
							xyznumdirs = Deconv_data_array[i]->GetODFMaxima(x,y,z)->GetNumberOfDirections();
							if (labels_exist) // which doesn't as of now
							{
							  if (0)
								  //(label_image_data->GetValue(x, y, z) > 0) //? The label_image_data has a value equal to one for a fanning configuration
								{
									//?
									xyznumdirs = DWI_data_array[0]->GetNSize();
								}
							}
	
							if (xyznumdirs > 3)
							{
								xyznumdirs = 3;
							}
	
							if ( xyznumdirs == 1)
							{
								dataInd = 0 + orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
	
								dir.x = Deconv_data_array[i]->GetODFMaxima(x,y,z)->GetXComponent(0);
								dir.y = Deconv_data_array[i]->GetODFMaxima(x,y,z)->GetYComponent(0);
								dir.z = Deconv_data_array[i]->GetODFMaxima(x,y,z)->GetZComponent(0);
	
								value = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);

								if (value < FLT_EPSILON)
								{
									mask_image_data->SetValue(x, y, z, 0);
								}
								else
								{
									dir.x = (float) dir.x/value;
									dir.y = (float) dir.y/value;
									dir.z = (float) dir.z/value;
		
									if (dataArray[dataInd] < 1)
									{
										dataArray[dataInd] = 1;
										ind = histogram_image_array[0]->FindClosestDirection(dir);
										offInd = orientNum*xsize*ysize*zsize;
										dataArray[dataInd + offInd*(ind+6)] += 1;
		
										dataArray[dataInd + offInd*(3)] = dir.x;
										dataArray[dataInd + offInd*(4)] = dir.y;
										dataArray[dataInd + offInd*(5)] = dir.z;
									}
									else
									{
										dataArray[dataInd] += 1;
										ind = histogram_image_array[0]->FindClosestDirection(dir);
										offInd = orientNum*xsize*ysize*zsize;
										dataArray[dataInd + offInd * (ind+6)] += 1;
		
										mdir[0].x = dataArray[dataInd + offInd*(3)];
										mdir[0].y = dataArray[dataInd + offInd*(4)];
										mdir[0].z = dataArray[dataInd + offInd*(5)];
		
										if ( dot(mdir[0], dir) < 0 )
										{
											dir.x = -1*dir.x;
											dir.y = -1*dir.y;
											dir.z = -1*dir.z;
										}
		
										mdir[0].x = mdir[0].x + (dir.x - mdir[0].x)/dataArray[dataInd];
										mdir[0].y = mdir[0].y + (dir.y - mdir[0].y)/dataArray[dataInd];
										mdir[0].z = mdir[0].z + (dir.z - mdir[0].z)/dataArray[dataInd];
		
										dataArray[dataInd + offInd*(3)] = mdir[0].x;
										dataArray[dataInd + offInd*(4)] = mdir[0].y;
										dataArray[dataInd + offInd*(5)] = mdir[0].z;
									}
								}								
							}
							else if (xyznumdirs == 2)
							{
								Dirs = new Directions(xyznumdirs);
								dataInd = orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
				
								if (dataArray[1 + dataInd] < 1)
								{
									for (j=0; j<xyznumdirs; j++)
									{
										dir = Deconv_data_array[i]->GetODFMaxima(x, y, z)->GetDirection(j);
										value = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
												
										dir.x = (float) dir.x/value;
										dir.y = (float) dir.y/value;
										dir.z = (float) dir.z/value;					
										Dirs->SetDirection(j, dir);
												
										dataArray[j + 1 + dataInd] = 1;
										ind = histogram_image_array[j + 1]->FindClosestDirection(Dirs->GetDirection(j));
										offInd = orientNum*xsize*ysize*zsize;
										dataArray[j + 1 + dataInd + offInd * (ind+6)] += 1;
				
										dataArray[j + 1 + dataInd + offInd*(3)] = Dirs->GetDirection(j).x;
										dataArray[j + 1 + dataInd + offInd*(4)] = Dirs->GetDirection(j).y;
										dataArray[j + 1 + dataInd + offInd*(5)] = Dirs->GetDirection(j).z;				
									}
								}
								else
								{
									for (j=0; j<xyznumdirs; j++)
									{
										dir = Deconv_data_array[i]->GetODFMaxima(x, y, z)->GetDirection(j);
										value = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
				
										dir.x = (float) dir.x/value;
										dir.y = (float) dir.y/value;
										dir.z = (float) dir.z/value;
										Dirs->SetDirection(j, dir);
				
										offInd = orientNum*xsize*ysize*zsize;
										mdir[j].x = dataArray[j + 1 + dataInd + offInd*(3)];
										mdir[j].y = dataArray[j + 1 + dataInd + offInd*(4)];
										mdir[j].z = dataArray[j + 1 + dataInd + offInd*(5)];
									}
				
									maxDot[0] = -1*FLT_MAX;
									for (j=0; j<xyznumdirs; j++)
									{
										for (k=0; k<xyznumdirs; k++)
										{
											dot_matrix[j][k] = float_abs(dot(Dirs->GetDirection(j), mdir[k]));
				
											if (maxDot[0] < dot_matrix[j][k])
											{
												maxDot[0] = dot_matrix[j][k];
												sInd[0] = j;
												mInd[0] = k;
											}
										}
									}
				
									sInd[1] = xyznumdirs - sInd[0] -1;
									mInd[1] = xyznumdirs - mInd[0] -1;
				
									for (j=0; j<xyznumdirs; j++)
									{
										dataArray[mInd[j] + 1 + dataInd] += 1;
										ind = histogram_image_array[mInd[j] + 1]->FindClosestDirection(Dirs->GetDirection(sInd[j]));
										offInd = orientNum*xsize*ysize*zsize;
										dataArray[mInd[j] + 1 + dataInd + offInd * (ind+6)] += 1;
				
										dir = Dirs->GetDirection(sInd[j]);
										if ( dot(mdir[mInd[j]], Dirs->GetDirection(sInd[j])) < 0 )
										{
											dir.x = -1*dir.x;
											dir.y = -1*dir.y;
											dir.z = -1*dir.z;
				
											Dirs->SetDirection(sInd[j], dir);
										}
				
										mdir[mInd[j]].x = mdir[mInd[j]].x + (Dirs->GetDirection(sInd[j]).x - mdir[mInd[j]].x)/dataArray[mInd[j] + 1 + dataInd];
										mdir[mInd[j]].y = mdir[mInd[j]].y + (Dirs->GetDirection(sInd[j]).y - mdir[mInd[j]].y)/dataArray[mInd[j] + 1 + dataInd];
										mdir[mInd[j]].z = mdir[mInd[j]].z + (Dirs->GetDirection(sInd[j]).z - mdir[mInd[j]].z)/dataArray[mInd[j] + 1 + dataInd];
			
										dataArray[mInd[j] + 1 + dataInd + offInd*(3)] = mdir[mInd[j]].x;
										dataArray[mInd[j] + 1 + dataInd + offInd*(4)] = mdir[mInd[j]].y;
										dataArray[mInd[j] + 1 + dataInd + offInd*(5)] = mdir[mInd[j]].z;
									}
								}
				
								Dirs->Release();
							}
							else if (xyznumdirs == 3)
							{
								Dirs = new Directions(xyznumdirs);
								dataInd = orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
	
								if (dataArray[3 + dataInd] < 1)
								{
									for (j=0; j<xyznumdirs; j++)
									{
										dir = Deconv_data_array[i]->GetODFMaxima(x, y, z)->GetDirection(j);
										value = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
	
										dir.x = (float) dir.x/value;
										dir.y = (float) dir.y/value;
										dir.z = (float) dir.z/value;
										Dirs->SetDirection(j, dir);
	
										dataArray[j + 3 + dataInd] = 1;
										ind = histogram_image_array[j + 3]->FindClosestDirection(Dirs->GetDirection(j));
										offInd = orientNum*xsize*ysize*zsize;
										dataArray[j + 3 + dataInd + offInd * (ind+6)] += 1;
	
										dataArray[j + 3 + dataInd + offInd*(3)] = Dirs->GetDirection(j).x;
										dataArray[j + 3 + dataInd + offInd*(4)] = Dirs->GetDirection(j).y;
										dataArray[j + 3 + dataInd + offInd*(5)] = Dirs->GetDirection(j).z;
									}
								}
								else
								{
									for (j=0; j<xyznumdirs; j++)
									{
										dir = Deconv_data_array[i]->GetODFMaxima(x, y, z)->GetDirection(j);
										value = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
	
										dir.x = (float) dir.x/value;
										dir.y = (float) dir.y/value;
										dir.z = (float) dir.z/value;
										Dirs->SetDirection(j, dir);
	
										offInd = orientNum*xsize*ysize*zsize;
										mdir[j].x = dataArray[j + 3 + dataInd + offInd*(3)];
										mdir[j].y = dataArray[j + 3 + dataInd + offInd*(4)];
										mdir[j].z = dataArray[j + 3 + dataInd + offInd*(5)];
									}
	
									maxDot[0] = maxDot[1] = -1*FLT_MAX;
	
									for (j=0; j<xyznumdirs; j++)
									{
										for (k=0; k<xyznumdirs; k++)
										{
											dot_matrix[j][k] = float_abs(dot(Dirs->GetDirection(j), mdir[k]));
	
											if (maxDot[0] < dot_matrix[j][k])
											{
												maxDot[0] = dot_matrix[j][k];
												sInd[0] = j;
												mInd[0] = k;
											}
										}
									}
	
									for (j=0; j<xyznumdirs; j++)
									{
										for (k=0; k<xyznumdirs; k++)
										{
											if ((j != sInd[0]) && (k != mInd[0]))
											{
												dot_matrix[j][k] = float_abs(dot(Dirs->GetDirection(j), mdir[k]));
	
												if (maxDot[1] < dot_matrix[j][k])
												{
													maxDot[1] = dot_matrix[j][k];
													sInd[1] = j;
													mInd[1] = k;
												}
	
											}
										}
									}
	
									sInd[2] = xyznumdirs - sInd[0] - sInd[1];
									mInd[2] = xyznumdirs - mInd[0] - mInd[1];
	
									for (j=0; j<xyznumdirs; j++)
									{
										dataArray[mInd[j] + 3 + dataInd] += 1;
										ind = histogram_image_array[mInd[j] + 3]->FindClosestDirection(Dirs->GetDirection(sInd[j]));
										offInd = orientNum*xsize*ysize*zsize;
										dataArray[mInd[j] + 3 + dataInd + offInd * (ind+6)] += 1;
	
										dir = Dirs->GetDirection(sInd[j]);
										if ( dot(mdir[mInd[j]], dir) < 0 )
										{
											dir.x = -1*dir.x;
											dir.y = -1*dir.y;
											dir.z = -1*dir.z;
	
											Dirs->SetDirection(sInd[j], dir);
										}
	
										mdir[mInd[j]].x = mdir[mInd[j]].x + (Dirs->GetDirection(sInd[j]).x - mdir[mInd[j]].x)/dataArray[mInd[j] + 3 + dataInd];
										mdir[mInd[j]].y = mdir[mInd[j]].y + (Dirs->GetDirection(sInd[j]).y - mdir[mInd[j]].y)/dataArray[mInd[j] + 3 + dataInd];
										mdir[mInd[j]].z = mdir[mInd[j]].z + (Dirs->GetDirection(sInd[j]).z - mdir[mInd[j]].z)/dataArray[mInd[j] + 3 + dataInd];
	
										dataArray[mInd[j] + 3 + dataInd + offInd*(3)] = mdir[mInd[j]].x;
										dataArray[mInd[j] + 3 + dataInd + offInd*(4)] = mdir[mInd[j]].y;
										dataArray[mInd[j] + 3 + dataInd + offInd*(5)] = mdir[mInd[j]].z;
									}
								}
	
								Dirs->Release();
							}
							//else if (xyznumdirs == DWI_data_array[0]->GetNSize()) //faning
							//{

							//}


						
						max_for_norm=Deconv_data_array[i]->GetMaximumValue(x,y,z);
						
						//cout << "set ODF values\n";
						for (n=0; n<DWI_data_array[0]->GetNSize(); n++)
						  {
						    //cout << Deconv_data_array[i]->GetValue(x,y,z,n) << " " << (iterations*max_for_norm) << endl;
						    full_avg_odf->SetValue(x,y,z,n,full_avg_odf->GetValue(x,y,z,n)+(Deconv_data_array[i]->GetValue(x,y,z,n)*value_for_vol_one/(iterations*max_for_norm)));
						  }
						}
	


					}
				}
			}
			//cout << i<< endl;
			Deconv_data_array[i]->Release();
		}
	
		//Normalise the number of times a label is observed and also normalise the mean direction to a unit vector.
		//Compute the histogram values for each fibre maximum by the number of time that maximum is observed and the normalization factor
		//which will be used next to normalize the histograms to unit volume
	

		if (!justpdf)
		  {



		Dirs = DWI_data_array[0]->GetDirections();
		for (i=0; i<orientNum; i++)
		{
			for (x=0; x<xsize; x++)
			{
				for (y=0; y<ysize; y++)
				{
					for (z=0; z<zsize; z++)
					{
						if (mask_image_data->GetValue(x,y,z) > maskVal)
						{
							dataInd = i + orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
	
							if (dataArray[dataInd] != 0)
							{
								volume[i]  = 0;
								for (j=0; j<Dirs->GetNumberOfDirections(); j++)
								{
									dataArray[dataInd + offInd*(j+6)] = (float)dataArray[dataInd + offInd*(j+6)]/(dataArray[dataInd]);
									volume[i] += (4 * M_PI * (float)pow(dataArray[dataInd + offInd*(j+6)], 3))/(3 * Dirs->GetNumberOfDirections());
								}
	
								offInd = orientNum*xsize*ysize*zsize;
								mdir[0].x = dataArray[dataInd + offInd*(3)];
								mdir[0].y = dataArray[dataInd + offInd*(4)];
								mdir[0].z = dataArray[dataInd + offInd*(5)];
	
								value = sqrt(pow(mdir[0].x, 2) + pow(mdir[0].y, 2) + pow(mdir[0].z, 2));
	
								if (value != 0)
								{
									mdir[0].x = (float) mdir[0].x/value;
									mdir[0].y = (float) mdir[0].y/value;
									mdir[0].z = (float) mdir[0].z/value;
								}
	
								if (isnan(mdir[0].x) || isnan(mdir[0].y) || isnan(mdir[0].z))
								{
									cout << "NAN!\n";
								}
	
								dataArray[dataInd + offInd*(3)] = mdir[0].x;
								dataArray[dataInd + offInd*(4)] = mdir[0].y;
								dataArray[dataInd + offInd*(5)] = mdir[0].z;
	
								value = (float) dataArray[dataInd]/iterations;
								dataArray[dataInd] = value;
	
								float val1 = pow(volume[i], (float)(1.0/3.0));
								for (j=0; j<Dirs->GetNumberOfDirections(); j++)
								{
									dataArray[dataInd + offInd*(j+6)] = (float)dataArray[dataInd + offInd*(j+6)]/val1;
								}
							}
						}
					}
				}
			}
		}
	
		// Now, the maxima of different fibre configurations are compared, and true confidence values will be assigned to them.
		//
		// M_PI*5/180: Maximim allowed value for angular difference between two given mean directions. (Hardcoded for now)
	
		offInd = orientNum*xsize*ysize*zsize;
	
		for (x=0; x<xsize; x++)
		{
			for (y=0; y<ysize; y++)
			{
				for (z=0; z<zsize; z++)
				{
					dataInd = orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
					if (mask_image_data->GetValue(x,y,z) > maskVal)
					{
						for (i=0; i<3; i++)
						{
							for (m=0; m<num_of_dir[i]; m++)
							{
								for (j=0; j<orientNum; j++)
								{
									labelsCheck[j] = false;
								}
	
								labelCount = dataArray[i + int_max(0, i-1) + m + dataInd];
								if (labelCount > 0)
								{
									if (dataArray[i + int_max(0, i-1) + m + dataInd + offInd *(1)] == 0)
									{
										mdir[i].x = dataArray[i + int_max(0, i-1) + m + dataInd + offInd *(3)];
										mdir[i].y = dataArray[i + int_max(0, i-1) + m + dataInd + offInd *(4)];
										mdir[i].z = dataArray[i + int_max(0, i-1) + m + dataInd + offInd *(5)];
	
										for (j=i+1; j<3; j++)
										{
											maxDot[0] = -1*FLT_MAX;
											ind = -1;
											for (n=0; n<num_of_dir[j]; n++)
											{
												labelCount = dataArray[j + int_max(0, j-1) + n + dataInd];
												if (labelCount > 0)
												{
													if (dataArray[j + int_max(0, j-1) + n + dataInd + offInd *(1)] == 0)
													{
														mdir[j].x = dataArray[j + int_max(0, j-1) + n + dataInd + offInd *(3)];
														mdir[j].y = dataArray[j + int_max(0, j-1) + n + dataInd + offInd *(4)];
														mdir[j].z = dataArray[j + int_max(0, j-1) + n + dataInd + offInd *(5)];
	
														dotProd = float_abs(dot(mdir[i], mdir[j]));
														value = sqrt(pow(mdir[i].x, 2) + pow(mdir[i].y, 2) + pow(mdir[i].z, 2));
														value = value * sqrt(pow(mdir[j].x, 2) + pow(mdir[j].y, 2) + pow(mdir[j].z, 2));
														dotProd = (float) dotProd/value;
	
														if (maxDot[0] < dotProd)
														{
															maxDot[0] = dotProd;
															ind = n;
														}
													}
												}
											}
	
											if (ind != -1)
											{
												if (maxDot[0] > cos((float)(M_PI*10)/180))
												{
													value = dataArray[j + int_max(0, j-1) + ind + dataInd];
													if (dataArray[i + int_max(0, i-1) + m + dataInd + offInd*(1)] == 0)
													{
														value += dataArray[i + int_max(0, i-1) + m + dataInd];
													}
													else
													{
														value += dataArray[i + int_max(0, i-1) + m + dataInd + offInd*(1)];
													}
													dataArray[i + int_max(0, i-1) + m + dataInd + offInd*(1)] = value;
	
													labelsCheck[j + int_max(0, j-1) + ind] = true;
												}
											}
										}
	
										//if (labels_exist) //fanning
										//{
										//}
	
										for (j=i+1; j<3; j++)
										{
											for (n=0; n<num_of_dir[j]; n++)
											{
												if (labelsCheck[j+int_max(0, j-1)+n] == true)
												{
													value = dataArray[i + int_max(0, i-1) + m + dataInd + offInd*(1)];
													dataArray[j + int_max(0, j-1) + n + dataInd + offInd*(1)]= value;
												}
											}
										}
									}
								}
	
								if (dataArray[i + int_max(0, i-1) + m + dataInd + offInd *(1)] == 0)
								{
									dataArray[i + int_max(0, i-1) + m + dataInd + offInd *(1)] = dataArray[i + int_max(0, i-1) + m + dataInd];
								}
							}
						}
	
						//if (labels_exist) //fanning
						//{
						//	dataArray[6 + dataInd + offInd *(1)] = dataArray[6 + dataInd];
						//}
					}
				}
			}
		}
	
		//Here the 95% confidence interval is computed of all except the fanning configuration:
	
		float upperPercent, alpha, cfInterval;
		float dirProd[Dirs->GetNumberOfDirections()];
		int indArray[Dirs->GetNumberOfDirections()];
	
		offInd = orientNum*xsize*ysize*zsize;
	
		for (x=0; x<xsize; x++)
		{
			for (y=0; y<ysize; y++)
			{
				for (z=0; z<zsize; z++)
				{
					if (mask_image_data->GetValue(x,y,z) > maskVal)
					{
						for (i=0; i<6; i++)
						{
							dataInd = i + orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
							if (dataArray[dataInd] > 0)
							{
								dir.x = dataArray[dataInd + offInd*(3)];
								dir.y = dataArray[dataInd + offInd*(4)];
								dir.z = dataArray[dataInd + offInd*(5)];
	
								ind = histogram_image_array[i]->FindClosestDirection(dir);
	
								for (j=0; j<Dirs->GetNumberOfDirections(); j++)
								{
									dir2 = Dirs->GetDirection(j);
									dirProd[j] = float_abs(dot(dir, dir2));
	
									value = sqrt(dir2.x*dir2.x + dir2.y*dir2.y + dir2.z*dir2.z);
									dirProd[j] = dirProd[j]/value;
									indArray[j] = j;
								}
	
								//sorting
								int max;
								for (j=0; j<Dirs->GetNumberOfDirections()-1; j++)
								{
									max = j;
									for (k=j+1; k<Dirs->GetNumberOfDirections(); k++)
									{
										if (dirProd[k] > dirProd[max])
										{
											max = k;
										}
									}
	
									value = indArray[j];
									indArray[j] = indArray[max];
									indArray[max] = (int)value;
	
									value = dirProd[j];
									dirProd[j] = dirProd[max];
									dirProd[max] = value;
								}
	
								upperPercent = 0;
								j = -1;
								// For the probability distribution of a random variable X, the θ percentage point (or lower percentage point) of the distribution
								// is x1, such that P(X<x1)=θ/100. The upper percentage point of the distribution is x2, such that P(X>x2)=θ/100.
								//while ((upperPercent < 0.95*dataArray[dataInd]*iterations)&& (j < Dirs->GetNumberOfDirections()-1))
							    while ((upperPercent < 0.95)&& (j < Dirs->GetNumberOfDirections()-1))
								{
									j++;
									upperPercent += (float)pow(dataArray[dataInd + offInd*(indArray[j]+6)], 3)* 4 * M_PI/(3 * Dirs->GetNumberOfDirections());
									// upperPercent += (float) dataArray[dataInd + offInd*(indArray[j]+6)];
								}
	
								//if (upperPercent < 0.95*dataArray[dataInd]*iterations)
								if (upperPercent < 0.95)
								{
									dataArray[dataInd + offInd*(2)] = (float)M_PI/2;
								}
								else
								{
									if (dirProd[j]>1)
									{
										dirProd[j] = 1;
									}
	
									alpha = acos(dirProd[j]);
									cfInterval = acos(1 - alpha/(2 * dataArray[dataInd] * iterations));
	
									if (isnan(cfInterval))
									{
									  //cout << "Why NAN?\n";
									}
	
									dataArray[dataInd + offInd*(2)] = cfInterval;
								}
							}
						}
					}
				}
			}
		}

		//Write the computed values into the final imagedata files
	
		for (x=0; x<xsize; x++)
		{
			for (y=0; y<ysize; y++)
			{
				for (z=0; z<zsize; z++)
				{
					if (mask_image_data->GetValue(x,y,z) > maskVal)
					{
						for (i=0; i<orientNum; i++)
						{
							ind = i + orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z;
							for (n=0; n<img_geom.n_size; n++)
							{
								dataInd = i + orientNum*x + orientNum*xsize*y + orientNum*xsize*ysize*z + orientNum*xsize*ysize*zsize*n;
	
								if (n<6)
								{
									value = dataArray[dataInd];
								}
								else
								{
									//value = (float) dataArray[dataInd]/(dataArray[ind] * iterations);
									value = dataArray[dataInd];
								}
	
								histogram_image_array[i]->SetValue(x, y, z, n, value);
							}
						}
					}
				}
			}
		}
		  }
	
		// Each histogram is next normalised to unit volume
	
		char inttostring[5];
		char OutFile[200];

		full_avg_odf->WriteMincFile("full_pdf.mnc",clobber, history,false);
	
		for (i=0; i<orientNum; i++)
		{
			memset(OutFile, '\0', 200);
			strcpy(OutFile, output_name);
			sprintf(inttostring,"%i",i);
			strcat(OutFile, inttostring);
			strcat(OutFile, ".mnc");
			histogram_image_array[i]->WriteMincFile(OutFile, clobber, history, false);
		}

		
	
		DWI_data_array[0]->Release();

		//response->Release();
		//response_voxel->Release();
	
		//if (remove_zero_b && algorithm!=FORECAST)
		//{
		//	response_nonzero_b->Release();
		//}
	
		//if (remove_zero_b)
		//{
		//	DWI_image_data->Release();
		//}
	
		SH_fit->Release();
		residuals->Release();
	
		free(Deconvolver_array);
		free(DWI_data_array);
		free(Deconv_data_array);
	
		//if (algorithm == FORECAST)
		//{
		//	lambda1->Release();
		//	lambda2->Release();
		//	}
	
		//if (use_e1)
		//{
		//	e1x->Release();
		//	e1y->Release();
		//	e1z->Release();
		//}
	
		mask_image_data->Release();
	
		for (i=0; i<6; i++)
		{
		  histogram_image_array[i]->Release();
		}
		free(histogram_image_array);
		free(dataArray);
		
	}
	
	void Usage()
	{
	  printf("Usage: mincProbDeconvolve [options] -DWIs <DWIs.mnc> -mask <calc_mask.mnc>  -o <out> \n\nDWIs: minc file must include a single nonzero bvalue series of DWIs in different directions and at least one b=0 image, which must be at the *beginning* of the series\nmask: voxels in which deconvolution will be performed\nout: base name (`no .mnc!') for the six output files, described below.\n\noptions:\n\n-response_voxel <MaskforFiberResponse.mnc>: a file in which the voxel or voxels that are considered the single fibre response function are nonzero.\n-iterations <number of iterations> default: 1000\n-lambda1 <lambda1.mnc> -lambda2 <lambda2.mnc>: eigenvalues (minc files) from a previous tensor fit: used for the FORECAST algorithm (The FORECAST algorithm is the only algorithm currently available - default.)\n-scalars_for_response <So> <l1> <l2>\n-clobber: overwrite existing <out.mnc>\n-order <Order of deconvolution> (default (and maximum) is 8)\n\n\nNOTE: all inputs must be sampled the same way\n\n\nOutput file format:\n\n");
	
	printf("The output consists of six 4D minc files.  The first contains the information for the single fibre geometry, the second and third files for the double crossing, and the fourth through sixth files contain the information for the triple crossing geometry.  For each of these files, the frames contain the following:\n\n");
	
	printf("\
	0: occurence of the geometry label (redundant across geometry)\n\
	1: occurence of the maximum in this file\n\
	2: 95%% confidence interval for the maximum (in radians)\n\
	3-5: mean direction\n\
	6-EOF:histogram of the orientations obtained with bootstrap procedure: orientations used for sampling are in the minc header.\n\
	\n");
	
	
	/*
	-sigma_theta: write out individual .txt files for the uncertainties associated with each vector in each geometry (one file per geometry): this input should be the base filename with no extension \n-maxima: write out individual .txt files for each vector direction in each geometry (one file per geometry): this input should be the base filename with no extension  \n
	*/
	
	exit(0);
	
	}
	
	int main(int argc, char *argv[])
	{
		int i,j;
		int iterations=1000;
	
		int order=8;
	
		char *input_name;
		char *mask_name;
		char *response_name;
		char *response_mask_name;
		char *output_name;
		char *sigma_file=NULL;
		char *maxima_file=NULL;
	
		bool clobber = false;
		bool use_e1 = false;
	
		char *e1x_name;
		char *e1y_name;
		char *e1z_name;
		char *lambda1_name=NULL;
		char *lambda2_name;
	
		char tempdirname[200];
		char *intempfilename_array[100];
		char *outtempfilename_array[100];
		char *masktempfilename_array[100];
		char *responsetempfilename_array[100];
		char *l1tempfilename_array[100];
		char *l2tempfilename_array[100];
	
		int slices_at_once=2;
		char command[1000];
		int zdim,ydim,tdim,xdim;
		float zstep;
	
		float l1,l2,So;
		So=0;

		if (argc < 3)//actually a lot more than 3 are required...
		{
			Usage();
		}
	
		// parse options: can be in any order
	
		for(i = 1; i < argc; i++)
		{
			if(!strcmp(argv[i], "-help"))
			{
				Usage();
			}
			else if(!strcmp(argv[i], "-clobber"))
			{
				clobber = true;
			}
			else if (!strcmp(argv[i], "-DWIs"))
			{
				if(++i >= argc) Usage();
				input_name = argv[i];
			}
			else if (!strcmp(argv[i], "-mask"))
			{
				if(++i >= argc) Usage();
				mask_name = argv[i];
			}
			/*
			else if (!strcmp(argv[i], "-response"))
			{
				if(++i >= argc) Usage();
				response_name = argv[i];
			}
			*/
			else if (!strcmp(argv[i], "-response_voxel"))
			{
				if(++i >= argc) Usage();
				response_mask_name = argv[i];
			}
			else if (!strcmp(argv[i], "-order"))
			{
				if(++i >= argc) Usage();
				order = atoi(argv[i]);
	
				if (order>8)
				{
					Error("Implemented only up to order 8 so far...");
				}
			}
			else if (!strcmp(argv[i], "-o"))
			{
				if(++i >= argc) Usage();
				output_name = argv[i];
			}
			else if (!strcmp(argv[i], "-iterations"))
			{
				if(++i >= argc) Usage();
				iterations = atoi(argv[i]);
			}
			/*
			else if (!strcmp(argv[i], "-sigma_theta"))
			{
				if(++i >= argc) Usage();
				sigma_file = argv[i];
			}
			else if (!strcmp(argv[i], "-maxima"))
			{
				if(++i >= argc) Usage();
				maxima_file = argv[i];
			}
			*/
			/*
			else if(!strcmp(argv[i], "-noremove_zero_b"))
			{
				remove_zero_b = false;
			}
	
			else if (!strcmp(argv[i], "-e1"))
			{
				if(++i >= argc) Usage();
				e1x_name = argv[i++];
				e1y_name = argv[i++];
				e1z_name = argv[i];
				use_e1=true;
			}
			*/
			else if (!strcmp(argv[i], "-lambda1"))
			{
				if(++i >= argc) Usage();
				lambda1_name = argv[i];
			}
			else if (!strcmp(argv[i], "-lambda2"))
			{
				if(++i >= argc) Usage();
				lambda2_name = argv[i];
			}
			else if (!strcmp(argv[i], "-scalars_for_response"))
			  {
			    if(++i >= argc) Usage();
			    So = (float) atof(argv[i++]);
			    l1 = (float) atof(argv[i++]);
			    l2 = (float) atof(argv[i]);
			  }
			else if(argv[i][0] == '-')
			{
				Usage();
			}
		}
	
		char *history=CreateHistoryString(argc, argv);
	
		
	
		//get the lambda1 and lambda2 constants from the input full files: make a dummy Deconvolver just to set and get these values (could reimplment set function here, but this eliminates repeated code):
	
	
	
		Deconvolver *deconvolver;
		ImageData *response;
		ImageData *l1im;
		ImageData *l2im;
		ImageData *resp;
		
	
		response=new ImageData(input_name);
		xdim=response->GetXSize();
		ydim=response->GetYSize();
		zdim=response->GetZSize();
		tdim=response->GetNSize();
		zstep=response->GetZStep();
	
		deconvolver=new Deconvolver();

		if (lambda1_name!=NULL)
		  {
		    
		    l1im=new ImageData(lambda1_name);
		    l2im=new ImageData(lambda2_name);
		    resp=new ImageData(response_mask_name);
		    deconvolver->SetResponseFunction(response);
		    deconvolver->SetLambdas(l1im,l2im,resp);
		    deconvolver->GetEigenvaluesAndSo(l1,l2,So);
		    delete(deconvolver);
		//these get deleted in deconvolver because we didn't retain them:
		//response->Release();
		//l1im->Release();
		//l2im->Release();
		//resp->Release();
		  }
		//none of the following deconvolver stuff is necessary: we just need the l1,l2,So values.
		else if (So>FLT_EPSILON)
		  {

		    deconvolver->SetLambdas(l1,l2,So);
		    delete(deconvolver);
		  }
		
		else//defaults: not a good idea!!  works for 2x2x2 on 3T acquired at MNI with our std protocol
		  {
		    //cout << "using defaults\n";
		    deconvolver->SetLambdasDefault();
		    l1=0.001814152052;
		    l2=0.0002350908912;
		    So=156.5972435;
		    delete(deconvolver);
	
		  }
	
	
		//do this in sections because of memory problems:
	
		bool break_up=false;
	
		if (break_up)
		  {
	
		sprintf(tempdirname,"/tmp/tmp%i",(int) time(NULL));
		sprintf(command,"mkdir %s\n",tempdirname);
		system(command);
	
		for (i=0;i<ceil(zdim/slices_at_once);i++)
		  {
		    j=slices_at_once*i;
		    intempfilename_array[i]=(char *) malloc(200*sizeof(char));
		    sprintf(intempfilename_array[i],"%s/tmpin%i.mnc",tempdirname,(int) time(NULL));
		    masktempfilename_array[i]=(char *) malloc(200*sizeof(char));
		    sprintf(masktempfilename_array[i],"%s/tmpmask%i.mnc",tempdirname,(int) time(NULL));
		    sprintf(command,"mincreshape -dimrange time=0,%i -dimrange zspace=%i,%i -dimrange yspace=0,%i -dimrange xspace=0,%i %s %s",tdim,j,slices_at_once,ydim,xdim,input_name,intempfilename_array[i]);
		    system(command);
		    sprintf(command,"mincreshape -dimrange zspace=%i,%i -dimrange yspace=0,%i -dimrange xspace=0,%i %s %s",j,slices_at_once,ydim,xdim,mask_name,masktempfilename_array[i]);
		    system(command);
		    outtempfilename_array[i]=(char *) malloc(200*sizeof(char));
		    sprintf(outtempfilename_array[i],"%s/tmpout%i.mnc",tempdirname,(int) time(NULL));
		    CalculateMaximaAndSigma(intempfilename_array[i],masktempfilename_array[i], order, outtempfilename_array[i], iterations, l1, l2, So, clobber, history);
	
		  }
	
		if (clobber)
		  {
		    if (zstep<0)
		      {
			sprintf(command,"mincconcat -concat_dimension zspace -descending %s/tmpout* %s -clobber",tempdirname,output_name);
		      }
		    else
		      {
			sprintf(command,"mincconcat -concat_dimension zspace %s/tmpout* %s -clobber",tempdirname,output_name);
		      }
		  }
		else
		  {
		    if (zstep<0)
		      {
		    sprintf(command,"mincconcat -concat_dimension zspace -descending %s/tmpout* %s",tempdirname,output_name);
		      }
		    else
		      {
			sprintf(command,"mincconcat -concat_dimension zspace %s/tmpout* %s",tempdirname,output_name);
		      }
		  }
		system(command);
		if (clobber)
		  {
		    if (zstep<0)
		      {
			sprintf(command,"mincconcat -concat_dimension zspace -descending %s/tmpout* %s -clobber",tempdirname,output_name);
		      }
		    else
		      {
			sprintf(command,"mincconcat -concat_dimension zspace %s/tmpout* %s -clobber",tempdirname,output_name);
		      }
		  }
		else
		  {
		    if (zstep<0)
		      {
		    sprintf(command,"mincconcat -concat_dimension zspace -descending %s/tmpout* %s",tempdirname,output_name);
		      }
		    else
		      {
			sprintf(command,"mincconcat -concat_dimension zspace %s/tmpout* %s",tempdirname,output_name);
		      }
		  }
		system(command);
		// clean up; DO: clean up the arrays
		sprintf(command,"rm -rf %s\n",tempdirname);
		//system(command);
		  }
		else
		  {
	
		    CalculateMaximaAndSigma(input_name,mask_name, order, output_name, iterations, l1, l2, So, clobber, history);
		  }
	
		free(history);
		return(0);
	}
	
	
