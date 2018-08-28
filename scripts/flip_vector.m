function flip_vector(xfile,yfile,zfile)
% 
%%% Ensure that all vectors within a neighborhood have minimum angle between them
%% (this comes about because diffusion data has no directionality in the sense that opposite directions cannot be resolved)

%%%%%%%%%%
%% Ilana 2007/04
%%%%%%%%%
%%
% function [x_flip,y_flip,z_flip] = flip_vector(xfile,yfile,zfile,maskf)
%    times, which, start)
%
% 	x: x component of vector (3D)
% 	y: y component of vector (3D)
% 	z: z component of vector (3D)
%%% Outputs
%	x_flip : x component either flipped or not
%	y_flip : y component either flipped or not
%	z_flip : z component either flipped or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementation:
%%	For each voxel compare its direction with the 26 (because 3D) around it. The dot product between the current one
%%	and all around should be positive (indicates they're in the same direction). If not need flip the vectors that 
%%	are not aligned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% April 19th 2007 IL: removed mask option, multliply with mask before processing
%%

%% open images
%h_m = openimage([maskf]);
h_x = openimage([xfile]);
h_y = openimage([yfile]);
h_z = openimage([zfile]);


%% get info on volumes
num_slices = getimageinfo(h_x,'NumSlices');
height = getimageinfo(h_x,'ImageHeight');
width  = getimageinfo(h_x,'ImageWidth');

%MASK=getallimages(maskf);
X=getallimages(xfile);
Y=getallimages(yfile);
Z=getallimages(zfile);

x=reshape(X, width, height,num_slices);
y=reshape(Y, width, height,num_slices);
z=reshape(Z, width, height,num_slices);
%mask=reshape(MASK, width, height,num_slices);

%%initialize (once in matlab, the slices are displayed as heigth*width (ex. 256*96)
x_flip=zeros(height*width,num_slices); % new x vector file
y_flip=zeros(height*width,num_slices); % new y vector file
z_flip=zeros(height*width,num_slices); % new z vector file

r=[];
c=[];
s=[];
%for k=1:num_slices 
%   [r_s,c_s] = find(mask(:,:,k) > 1-0.5 & mask(:,:,k) < 1+0.5); %% get all voxels that aren't 0 
%   r=[r,r_s'];
%   c=[c,c_s'];
%   s=[s,k*ones(1,length(r_s))];
%end

for k=1:num_slices 
   [r_s,c_s] = find(x(:,:,k) > 1-0.5 & x(:,:,k) < 1+0.5); %% get all voxels that aren't 0 
   r=[r,r_s'];
   c=[c,c_s'];
   s=[s,k*ones(1,length(r_s))];
end


for i=1:length(c)
   %% have to check the 26 adjacent voxels (whether they're in the mask or not)
   %% flip if dot product is not positive
   
   %%%%%%%% if previous slice present
   if(s(i)-1 > 0)% previous slice
    if(dot([x(r(i),c(i),s(i)-1),y(r(i),c(i),s(i)-1),z(r(i),c(i),s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0)

       x(r(i),c(i),s(i)-1)= - x(r(i),c(i),s(i)-1);
       y(r(i),c(i),s(i)-1)= - y(r(i),c(i),s(i)-1);
       z(r(i),c(i),s(i)-1)= - z(r(i),c(i),s(i)-1);
    end %same, previous slice

    if(c(i)-1 > 0) % previous column   
      if(dot([x(r(i),c(i)-1,s(i)-1),y(r(i),c(i)-1,s(i)-1),z(r(i),c(i)-1,s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice, previous column

         x(r(i),c(i)-1,s(i)-1)= - x(r(i),c(i)-1,s(i)-1);
         y(r(i),c(i)-1,s(i)-1)= - y(r(i),c(i)-1,s(i)-1);
         z(r(i),c(i)-1,s(i)-1)= - z(r(i),c(i)-1,s(i)-1);
      end %% previous slice, previous column

      if(r(i)-1 > 0) 
        if(dot([x(r(i)-1,c(i)-1,s(i)-1),y(r(i)-1,c(i)-1,s(i)-1),z(r(i)-1,c(i)-1,s(i)-1)], [x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,previous col,previous row

           x(r(i)-1,c(i)-1,s(i)-1)= - x(r(i)-1,c(i)-1,s(i)-1);
           y(r(i)-1,c(i)-1,s(i)-1)= - y(r(i)-1,c(i)-1,s(i)-1);
           z(r(i)-1,c(i)-1,s(i)-1)= - z(r(i)-1,c(i)-1,s(i)-1);
        end % previous slice, previous column, previous row
      end

      if(r(i)+1 > 0) 
        if(dot([x(r(i)+1,c(i)-1,s(i)-1),y(r(i)+1,c(i)-1,s(i)-1),z(r(i)+1,c(i)-1,s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,previous col,next row

           x(r(i)+1,c(i)-1,s(i)-1)= - x(r(i)+1,c(i)-1,s(i)-1);
           y(r(i)+1,c(i)-1,s(i)-1)= - y(r(i)+1,c(i)-1,s(i)-1);
           z(r(i)+1,c(i)-1,s(i)-1)= - z(r(i)+1,c(i)-1,s(i)-1);
        end % previous slice, previous column, next row
      end
    end %previous column

    if(c(i)+1 > 0) % next column   
      if(dot([x(r(i),c(i)+1,s(i)-1),y(r(i),c(i)+1,s(i)-1),z(r(i),c(i)+1,s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice, next column

         x(r(i),c(i)+1,s(i)-1)= - x(r(i),c(i)+1,s(i)-1);
         y(r(i),c(i)+1,s(i)-1)= - y(r(i),c(i)+1,s(i)-1);
         z(r(i),c(i)+1,s(i)-1)= - z(r(i),c(i)+1,s(i)-1);
      end %% previous slice, next column

      if(r(i)-1 > 0) 
        if(dot([x(r(i)-1,c(i)+1,s(i)-1),y(r(i)-1,c(i)+1,s(i)-1),z(r(i)-1,c(i)+1,s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,next col,previous row

           x(r(i)-1,c(i)+1,s(i)-1)= - x(r(i)-1,c(i)+1,s(i)-1);
           y(r(i)-1,c(i)+1,s(i)-1)= - y(r(i)-1,c(i)+1,s(i)-1);
           z(r(i)-1,c(i)+1,s(i)-1)= - z(r(i)-1,c(i)+1,s(i)-1);
        end % previous slice, next column, previous row
      end

      if(r(i)+1 > 0) 
        if(dot([x(r(i)+1,c(i)+1,s(i)-1),y(r(i)+1,c(i)+1,s(i)-1),z(r(i)+1,c(i)+1,s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,next col,next row

           x(r(i)+1,c(i)+1,s(i)-1)= - x(r(i)+1,c(i)+1,s(i)-1);
           y(r(i)+1,c(i)+1,s(i)-1)= - y(r(i)+1,c(i)+1,s(i)-1);
           z(r(i)+1,c(i)+1,s(i)-1)= - z(r(i)+1,c(i)+1,s(i)-1);
        end % previous slice, next column, next row
      end

    end %next column

    if(r(i)-1 > 0) %previous row, current column
      if(dot([x(r(i)-1,c(i),s(i)-1),y(r(i)-1,c(i),s(i)-1),z(r(i)-1,c(i),s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,current col,previous row

         x(r(i)-1,c(i),s(i)-1)= - x(r(i)-1,c(i),s(i)-1);
         y(r(i)-1,c(i),s(i)-1)= - y(r(i)-1,c(i),s(i)-1);
         z(r(i)-1,c(i),s(i)-1)= - z(r(i)-1,c(i),s(i)-1);
      end % previous slice, current column, previous row
    end


    if(r(i)+1 > 0) %next row, current column
      if(dot([x(r(i)+1,c(i),s(i)-1),y(r(i)+1,c(i),s(i)-1),z(r(i)+1,c(i),s(i)-1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,current col,next row

         x(r(i)+1,c(i),s(i)-1)= - x(r(i)+1,c(i),s(i)-1);
         y(r(i)+1,c(i),s(i)-1)= - y(r(i)+1,c(i),s(i)-1);
         z(r(i)+1,c(i),s(i)-1)= - z(r(i)+1,c(i),s(i)-1);
      end % previous slice, current column, next row
    end
				
   end %%%%%%%%%%%%%%previous slice

%%%%%%%% if next slice present
   if(s(i)+1 > 0)% next slice
    if(dot([x(r(i),c(i),s(i)+1),y(r(i),c(i),s(i)+1),z(r(i),c(i),s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0)

       x(r(i),c(i),s(i)+1)= - x(r(i),c(i),s(i)+1);
       y(r(i),c(i),s(i)+1)= - y(r(i),c(i),s(i)+1);
       z(r(i),c(i),s(i)+1)= - z(r(i),c(i),s(i)+1);
    end %same, next slice

    if(c(i)-1 > 0) % previous column   
      if(dot([x(r(i),c(i)-1,s(i)+1),y(r(i),c(i)-1,s(i)+1),z(r(i),c(i)-1,s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next slice, previous column

         x(r(i),c(i)-1,s(i)+1)= - x(r(i),c(i)-1,s(i)+1);
         y(r(i),c(i)-1,s(i)+1)= - y(r(i),c(i)-1,s(i)+1);
         z(r(i),c(i)-1,s(i)+1)= - z(r(i),c(i)-1,s(i)+1);
      end %% next slice, previous column

      if(r(i)-1 > 0) 
        if(dot([x(r(i)-1,c(i)-1,s(i)+1),y(r(i)-1,c(i)-1,s(i)+1),z(r(i)-1,c(i)-1,s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next slice,previous col,previous row

           x(r(i)-1,c(i)-1,s(i)+1)= - x(r(i)-1,c(i)-1,s(i)+1);
           y(r(i)-1,c(i)-1,s(i)+1)= - y(r(i)-1,c(i)-1,s(i)+1);
           z(r(i)-1,c(i)-1,s(i)+1)= - z(r(i)-1,c(i)-1,s(i)+1);
        end % next slice, previous column, previous row
      end

      if(r(i)+1 > 0) 
        if(dot([x(r(i)+1,c(i)-1,s(i)+1),y(r(i)+1,c(i)-1,s(i)+1),z(r(i)+1,c(i)-1,s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next slice,previous col,next row

           x(r(i)+1,c(i)-1,s(i)+1)= - x(r(i)+1,c(i)-1,s(i)+1);
           y(r(i)+1,c(i)-1,s(i)+1)= - y(r(i)+1,c(i)-1,s(i)+1);
           z(r(i)+1,c(i)-1,s(i)+1)= - z(r(i)+1,c(i)-1,s(i)+1);
        end % next slice, previous column, next row
      end

    end %previous column

    if(c(i)+1 > 0) % next column   
      if(dot([x(r(i),c(i)+1,s(i)+1),y(r(i),c(i)+1,s(i)+1),z(r(i),c(i)+1,s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next slice, next column

         x(r(i),c(i)+1,s(i)+1)= - x(r(i),c(i)+1,s(i)+1);
         y(r(i),c(i)+1,s(i)+1)= - y(r(i),c(i)+1,s(i)+1);
         z(r(i),c(i)+1,s(i)+1)= - z(r(i),c(i)+1,s(i)+1);
      end %% next slice, next column

      if(r(i)-1 > 0) 
        if(dot([x(r(i)-1,c(i)+1,s(i)+1),y(r(i)-1,c(i)+1,s(i)+1),z(r(i)-1,c(i)+1,s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next slice,next col,previous row

           x(r(i)-1,c(i)+1,s(i)+1)= - x(r(i)-1,c(i)+1,s(i)+1);
           y(r(i)-1,c(i)+1,s(i)+1)= - y(r(i)-1,c(i)+1,s(i)+1);
           z(r(i)-1,c(i)+1,s(i)+1)= - z(r(i)-1,c(i)+1,s(i)+1);
        end % next slice, next column, previous row
      end

      if(r(i)+1 > 0) 
        if(dot([x(r(i)+1,c(i)+1,s(i)+1),y(r(i)+1,c(i)+1,s(i)+1),z(r(i)+1,c(i)+1,s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next slice,next col,next row

           x(r(i)+1,c(i)+1,s(i)+1)= - x(r(i)+1,c(i)+1,s(i)+1);
           y(r(i)+1,c(i)+1,s(i)+1)= - y(r(i)+1,c(i)+1,s(i)+1);
           z(r(i)+1,c(i)+1,s(i)+1)= - z(r(i)+1,c(i)+1,s(i)+1);
        end % next slice, next column, next row
      end

    end %previous column

    if(r(i)-1 > 0) %previous row, current column
      if(dot([x(r(i)-1,c(i),s(i)+1),y(r(i)-1,c(i),s(i)+1),z(r(i)-1,c(i),s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,current col,previous row

         x(r(i)-1,c(i),s(i)+1)= - x(r(i)-1,c(i),s(i)+1);
         y(r(i)-1,c(i),s(i)+1)= - y(r(i)-1,c(i),s(i)+1);
         z(r(i)-1,c(i),s(i)+1)= - z(r(i)-1,c(i),s(i)+1);
      end % previous slice, current column, previous row
    end


    if(r(i)+1 > 0) %next row, current column
      if(dot([x(r(i)+1,c(i),s(i)+1),y(r(i)+1,c(i),s(i)+1),z(r(i)+1,c(i),s(i)+1)],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous slice,current col,next row

         x(r(i)+1,c(i),s(i)+1)= - x(r(i)+1,c(i),s(i)+1);
         y(r(i)+1,c(i),s(i)+1)= - y(r(i)+1,c(i),s(i)+1);
         z(r(i)+1,c(i),s(i)+1)= - z(r(i)+1,c(i),s(i)+1);
      end % previous slice, current column, next row
    end
				
   end %%%%%%%%%%%%%%%%%%%%%%%%%% next slice

%%%%%%%% if previous column present
   if(c(i)-1 > 0)%
    if(dot([x(r(i),c(i)-1,s(i)),y(r(i),c(i)-1,s(i)),z(r(i),c(i),s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous col, same row

       x(r(i),c(i)-1,s(i))= - x(r(i),c(i)-1,s(i));
       y(r(i),c(i)-1,s(i))= - y(r(i),c(i)-1,s(i));
       z(r(i),c(i)-1,s(i))= - z(r(i),c(i)-1,s(i));
    end %previous col, same row

    if(r(i)-1 > 0) % previous row   
      if(dot([x(r(i)-1,c(i)-1,s(i)),y(r(i)-1,c(i)-1,s(i)),z(r(i)-1,c(i)-1,s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %%previous column, previous row

         x(r(i)-1,c(i)-1,s(i))= - x(r(i)-1,c(i)-1,s(i));
         y(r(i)-1,c(i)-1,s(i))= - y(r(i)-1,c(i)-1,s(i));
         z(r(i)-1,c(i)-1,s(i))= - z(r(i)-1,c(i)-1,s(i));
      end %previous column, previous row
    end

    if(r(i)+1 > 0) % next row 
      if(dot([x(r(i)+1,c(i)-1,s(i)),y(r(i)+1,c(i)-1,s(i)),z(r(i)+1,c(i)-1,s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% previous col,next row

         x(r(i)+1,c(i)-1,s(i))= - x(r(i)+1,c(i)-1,s(i));
         y(r(i)+1,c(i)-1,s(i))= - y(r(i)+1,c(i)-1,s(i));
         z(r(i)+1,c(i)-1,s(i))= - z(r(i)+1,c(i)-1,s(i));
      end % previous col,next row
    end

   end %%%%%%%%%%%%%%previous column

%%%%%%%% if next column present
   if(c(i)+1 > 0)%
    if(dot([x(r(i),c(i)+1,s(i)),y(r(i),c(i)+1,s(i)),z(r(i),c(i),s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next col, same row

       x(r(i),c(i)+1,s(i))= - x(r(i),c(i)+1,s(i));
       y(r(i),c(i)+1,s(i))= - y(r(i),c(i)+1,s(i));
       z(r(i),c(i)+1,s(i))= - z(r(i),c(i)+1,s(i));
    end %next col, same row

    if(r(i)-1 > 0) % previous row   
      if(dot([x(r(i)-1,c(i)+1,s(i)),y(r(i)-1,c(i)+1,s(i)),z(r(i)-1,c(i)+1,s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %%next column, previous row

         x(r(i)-1,c(i)+1,s(i))= - x(r(i)-1,c(i)+1,s(i));
         y(r(i)-1,c(i)+1,s(i))= - y(r(i)-1,c(i)+1,s(i));
         z(r(i)-1,c(i)+1,s(i))= - z(r(i)-1,c(i)+1,s(i));
      end %next column, previous row
    end

      if(r(i)+1 > 0) % next row 
        if(dot([x(r(i)+1,c(i)+1,s(i)),y(r(i)+1,c(i)+1,s(i)),z(r(i)+1,c(i)+1,s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% next col,next row

           x(r(i)+1,c(i)+1,s(i))= - x(r(i)+1,c(i)+1,s(i));
           y(r(i)+1,c(i)+1,s(i))= - y(r(i)+1,c(i)+1,s(i));
           z(r(i)+1,c(i)+1,s(i))= - z(r(i)+1,c(i)+1,s(i));
        end % next col,next row
      end

    end %%%%%%%%%%%%%%%%next column

%%%%%%%% if previous row present
    if(r(i)-1 > 0) % previous row   
      if(dot([x(r(i)-1,c(i),s(i)),y(r(i)-1,c(i),s(i)),z(r(i)-1,c(i),s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %%current column, previous row

         x(r(i)-1,c(i),s(i))= - x(r(i)-1,c(i),s(i));
         y(r(i)-1,c(i),s(i))= - y(r(i)-1,c(i),s(i));
         z(r(i)-1,c(i),s(i))= - z(r(i)-1,c(i),s(i));
      end %current column, previous row
    end %%%%%%%% previous row
%%%%%%%% if next row present
    if(r(i)+1 > 0) % next row 
      if(dot([x(r(i)+1,c(i),s(i)),y(r(i)+1,c(i),s(i)),z(r(i)+1,c(i),s(i))],[x(r(i),c(i),s(i)),y(r(i),c(i),s(i)),z(r(i),c(i),s(i))]) < 0) %% current col,next row

         x(r(i)+1,c(i),s(i))= - x(r(i)+1,c(i),s(i));
         y(r(i)+1,c(i),s(i))= - y(r(i)+1,c(i),s(i));
         z(r(i)+1,c(i),s(i))= - z(r(i)+1,c(i),s(i));
      end % current col,next row
   end %%%%%%%% next row


end %voxels


%keyboard
new_x=newimage('flip_x.mnc',[0 num_slices],[xfile]);
X=reshape(x,width*height,num_slices);
putimages(new_x,X, [1:num_slices]);

new_y=newimage('flip_y.mnc',[0 num_slices],[yfile]);
Y=reshape(y,width*height,num_slices);
putimages(new_y,Y, [1:num_slices]);

new_z=newimage('flip_z.mnc',[0 num_slices],[zfile]);
Z=reshape(z,width*height,num_slices);
putimages(new_z,Z, [1:num_slices]);


return
