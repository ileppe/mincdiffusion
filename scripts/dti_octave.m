%function DTI(BASE_IMAGES, [SIGMA])  
%
%
%Calculates diffusion tensor, eigenvectors, eigenvalues, and anisotropy
%index, from a series of base diffusion-weighted images. Saves (as minc
%function DTI(BASE_IMAGES, [SIGMA])  
%
%
%Calculates diffusion tensor, eigenvectors, eigenvalues, and anisotropy
%index, from a series of base diffusion-weighted images. Saves (as minc
%files) the eigenvectors, eigenvalues, trace, AI, (and, optionally,
%chi-square map for the diffusion tensor fit). Also saves
%color-coded (rgb) map of x,y,z components of principal
%eigenvector.  All vector outputs are in world space.
%
%BASE_IMAGES is a minc file containing the base images
%for the study, with images corresponding to different b values in the time
%dimension. 
%
%Conventions used are:
%
%Z increases from patient inferior to superior
%Y increases from patient posterior to anterior
%X increases from patient left to right
%
%Internal function bmatrix.m must be modified if b values change 
%
%This code has only been implemented and tested for transverse base images,
%
%Currently does only single-exponential fit; can be modified to save other 
%parameters
%
%
%SIGMA=mean intensity in noise * sqrt(2/pi) 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %% Modified Dec 2006; Ilana Leppert
%% -change AI o FA % -change trace to MD (mean diffusivity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hh = dti_octave(base_images, sigma)


more off;
warning off;

chisquare=1;

if (nargin<2)
  chisquare=0;
end


%preprocess 1st: motion correction difficult
%with DWIs, as different contrast per direction.  Should implement something 
% like Rohde, 2004, for robust eddy current corection and motion correction

%also do any correction for scanner drift; not sure how to robustly do this: 
%not doing it now.

%matlab initialization:

more off;
%warning off

THRESH=0;                             %brain (b=0) 100; phantom 300: need to
                                      %empirically find a good threshold
                                      %for each new acquisition protocol,

%% IL -change AI o FA
%% IL -change trace to MD (mean diffusivity)
if (chisquare)
  
  outputs=char({'MD','e1x','e1y','e1z','FA','e2x','e2y','e2z', 'e3x','e3y','e3z','lambda1','lambda2','lambda3','chi2'});
else
  outputs=char({'MD','e1x','e1y','e1z','FA','e2x','e2y','e2z','e3x','e3y','e3z','lambda1','lambda2','lambda3'});
end
%DO: make these completely user-selectable


 
%unzip input file if not done: (this isn't necessary; openimage will unzip it)
%unix(['gunzip -q ',base_images]);


%open minc images file: 
%%% use niak rather than emma tools
  
%f = openimage(base_images,'r');
[hdr,VOL]=niak_read_vol(base_images);
%slices=getimageinfo(f,'NumSlices');
slices=hdr.info.dimensions(3);
%width=getimageinfo(f,'ImageWidth');  
width=hdr.info.dimensions(1);
%height=getimageinfo(f,'ImageHeight');
height=hdr.info.dimensions(2);
%numimages=getimageinfo(f,'NumFrames');
numimages=hdr.info.dimensions(4);



[b,number_images,zeroimage]=bmatrixfromheader(hdr); 

%alter format of b for single-exponential calculation:

b=[b(:,1) b(:,2) b(:,3) 2*b(:,4) 2*b(:,5) 2*b(:,6) ones(number_images,1)];  


%DO: possibly get rid of this and mincresample everything to transverse:

%the following determines whether the slices are transverse, coronal, or
%sagittal (assuming there aren't HEIGHT slices!): DO: change it so more
%robust.  note: the rest of the code currently only works for transverse

if(hdr.info.dimension_order=='xyzt')
  orientation = 't';
  slicedir = 'zspace';
end

if ~(strcmp(slicedir,'zspace'))
  error('images must be transverse');
end

%% reshape back to a 2D matrix, easier for now (DO: change to make
%operations on volume (3D) data)

vol=reshape(VOL, width*height,slices,number_images);
numvoxels = height*width;
info=zeros(numvoxels,size(outputs,1),slices); %matrix to hold all calculated
%parameters; change size if you add parameters to outputs, eg ADC (x,y,z),
%e3 (x,y,z), D1, D2, and compartment fractions

for s=1:slices

  %% IL change getallimages to emma version getimages 
  %images=getimages(f,s,1:numimages);
  % don't use emma
  images=vol(:,s,:);
	
  if number_images~=numimages
    error('number of b matrices and number of frames in input image must be equal');
  end

  %eliminate noise before fit by thresholding b=0 image: comment the following
  %if you want information on noise in tensor and AI maps 



  %% Also protect against signal values identical to 0, causes the log operation to return 0, and no regression can be calculated
  %% set 0 voxel to highest precision close to 0

  noise=images<=THRESH;
  index=find(images(:,zeroimage)>0);
  images(noise)=1;

  
  %% Initialize some variables to speed up loop
  chisqu=zeros(numvoxels,1); 
  z = zeros(size(b,1),length(index));
  diff = zeros(6,numvoxels);

  %% take ln of signal
  %% octave doesn't like this, if 'index' is empty, it creates an empty 'z'
  %% matlab alos creates en empty array, but weirdly it tsays size 0xsize(b,1)
  if(length(index)==0)
    index=1; %% create a dummy, will compute in 1 vox, which is not ideal
  end
  z=-log(images(index,:));
  
  %note D = (b\z) is faster and equivalent to D=regress(z,b)
  %(it solves the overdetermined system in a least-squares sense); 
  %use regress for stats:
  %[D,BINT,R,RINT,stats] = regress(z,b,0.05);
  D = b\z'; 
  

   if (chisquare)
  %if (chisquare & size(D,2)>0) 	
    for i = 1:length(index)
      SSE=0;
          
      for j=1:size(b,1)
          factor=0;
          for k=1:6
              factor=factor-b((j),k)*D(k,i);        
          end
          

          %nonlinear:
          %
          factor=exp(factor);
          SSE=SSE+(images(index(i),j)-images(index(i),1)*factor)^2;
      end
      info(index(i),15,s)=SSE/sigma^2;
    end %index loop
  end %chisquare
  
  diff(:,index)=D(1:6,:); %D is in mm^2/s
  
  %DO: should save F and p value from stats as well.  


  %clear images;


  %for each voxel, construct D, calculate trace, find eigenvalues & vectors, find largest, calculate volume ratio. 
  %% can't vectorize the 'det' and 'eig', they only accept 2D arrays

  Dsquare=zeros(3);

  %for i=index'  causes errors, with or without the '.
  for i = 1:length(index)
  
    Dsquare(1,1)=diff(1,index(i));
    Dsquare(1,2)=diff(4,index(i));
    Dsquare(1,3)=diff(5,index(i));
    Dsquare(2,1)=diff(4,index(i));
    Dsquare(2,2)=diff(2,index(i));
    Dsquare(2,3)=diff(6,index(i));
    Dsquare(3,1)=diff(5,index(i));
    Dsquare(3,2)=diff(6,index(i));
    Dsquare(3,3)=diff(3,index(i));  
    
    info(index(i),1,s)=1E3/3*(Dsquare(1,1) + Dsquare(2,2) + Dsquare(3,3)); %trace 
    %will be in 1E-6 mm^2/ms
    if det(Dsquare)==0
	 info(index(i),:,s)=zeros(1,size(outputs,1));%D is zero matrix: no diffusion
    elseif isfinite(det(Dsquare)) 
    %DO: should first check for positive-definite D.
       [V,E]=eig(Dsquare);
       %DO: check for negatives here!! reject voxels with neg
       %(don't get any in regions with signal)
       %if (sum(abs([E(1,1) E(2,2) E(3,3)])~=[E(1,1) E(2,2) E(3,3)])>0)
       %  disp('negative eigenvalues');
       %end
       %also check for orthogonality:
       %sum(V(:,1).*V(:,2));
       %sum(V(:,2).*V(:,3));
       %sum(V(:,1).*V(:,3));
       %e=abs([E(1,1) E(2,2) E(3,3)]);
       e=[E(1,1) E(2,2) E(3,3)];
       [Y,I]=max(e);
       [Y,J]=min(e);
  	  if I==1
  	     e1=V(:,1); %principal direction of diffusion
  	     lambda1=e(1);
  	     if e(2)>e(3)
  		     lambda2=e(2);
  		     lambda3=e(3);
  		     e2=V(:,2);
  		     e3=V(:,3);
  	     else
  		     lambda2=e(3);
  		     lambda3=e(2);
  		     e2=V(:,3);
  		     e3=V(:,2);
	     end;
  	   end;
	
	  if I==2
	    e1=V(:,2); %principal direction of diffusion
	    lambda1=e(2);
	    if e(1)>e(3)
	  	    lambda2=e(1);
	  	    lambda3=e(3);
	  	    e2=V(:,1);
	  	    e3=V(:,3);
	    else
	  	    lambda2=e(3);
	  	    lambda3=e(1);
	  	    e2=V(:,3);
	  	    e3=V(:,1);
	    end;
	  end;
	  if I==3
	  	  e1=V(:,3); %principal direction of diffusion
	  	  lambda1=e(3);
	  	  if e(2)>e(1)
	  		  lambda2=e(2);
	  		  lambda3=e(1);
	  		  e2=V(:,2);
	  		  e3=V(:,1);
	  	  else
	  		  lambda2=e(1);
	  		  lambda3=e(2);
	  		  e2=V(:,1);
	  		  e3=V(:,2);
	  	  end;
	
	  end;
 
	%eigenvectors:
		
	

	
	info(index(i),2,s)=e1(1);
	info(index(i),3,s)=e1(2);
	info(index(i),4,s)=e1(3);
	
	
	info(index(i),6,s)=e2(1);
	info(index(i),7,s)=e2(2);
	info(index(i),8,s)=e2(3);
	
	info(index(i),9,s)=e3(1);
	info(index(i),10,s)=e3(2);
	info(index(i),11,s)=e3(3);
	
	%different AIs: uncomment the one you want
	
	%AI=1-volume ratio
        %info(i,5)=1-(lambda1*lambda2*lambda3/((lambda1+lambda2+lambda3)/3)^3);
        
	%AI=RA, from AL Alexander
	%lambda_bar=mean([lambda1 lambda2 lambda3]);
	%info(i,5)=sqrt((lambda1-lambda_bar)^2+(lambda2-lambda_bar)^2+(lambda3-lambda_bar)^2)/(sqrt(6)*lambda_bar);
	
	%AI=FA
	lambda_bar=mean([lambda1 lambda2 lambda3]);
	tmp=sqrt(3/2)*sqrt((lambda1-lambda_bar)^2+(lambda2-lambda_bar)^2+(lambda3-lambda_bar)^2)/sqrt(lambda1^2+lambda2^2+lambda3^2);
	%clamp to maximum of 1
        if(tmp>1)
          info(index(i),5,s)=1;
        else info(index(i),5,s)=tmp; end
    
	%Lattice Index (A_DD):
	%compares the anisotropic parts of the diffusion tensor in two neighbors 
	
	%LI(r_i)=sqrt(3/8)*sqrt(d(r_i):d(r)/D(r_i)/D(r))+3/4*d(r_i):d(r)/(sqrt(D(r):D(r))*sqrt(D(r_i):D(r_i));
	%LI(r)=sum(26 nbhd)(1/eucldist_i)*LI(r_i)/sum(i=1:K??)(1/eucldist_i)
	
	%d=deviatoric=D-trace/3*I.
	
	%eigenvalues:
	
	info(index(i),12,s)=lambda1;
	info(index(i),13,s)=lambda2;
	info(index(i),14,s)=lambda3;
	


     else
	   info(index(i),:,s)=zeros(1,size(outputs,1));	%D is inf
     end;
   end;


%keyboard
%clear info;
   
    
end %slices loop


%get name for output images:

name=base_images(1:size(base_images,2)-4);


%write minc images:
%note: enable output selection here..don't go by k

for k=1:size(outputs,1)
  %% octave inserts some unwanted whitespaces
  tmp=strtrim(outputs(k,:));
  name2=strcat(name,'_',tmp,'.mnc');
  if exist(name2) unix(['rm -f ',name2]); end;  %yikes: DO:should have error message!!!	  
  %h = newimage2(name2, [0 slices height width], base_images);
  %putimages(h,info(:,k,:),[1:slices]);
  %% use niak toolkit
  hdr.file_name=name2;
  vol_out=reshape(info(:,k,:),width,height,slices);
  niak_write_vol(hdr,vol_out);
end  



%% IL -change AI to FA
%% do this in minctensor.pl
%rgb_name=strcat(name,'_rgb.mnc');
%if exist(rgb_name) unix(['rm -f ',rgb_name]); end;  %yikes: DO:should have error message!!!	
%unix(['rgbvector ',name,'_e1x.mnc ',name,'_e1y.mnc ',name,'_e1z.mnc ',name,'_FA.mnc ',rgb_name]);

%keyboard

clear;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%internal functions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b,number_images,zeroimage]=bmatrixfromheader(hdr)

%we assume we want to use all of the images: numpoints>1 means we don't use the
%(assumes there is only one) b=0 image.  the b=0 image is zeroimage and it
%is used for checking mask, thresholding...

num = find(niak_find_str_cell(hdr.details.acquisition.varatts,'bvalues'));
bvalues = hdr.details.acquisition.attvalue{num};

%these are minc x,y,z
num = find(niak_find_str_cell(hdr.details.acquisition.varatts,'direction_x'));
directionx = hdr.details.acquisition.attvalue{num};

num = find(niak_find_str_cell(hdr.details.acquisition.varatts,'direction_y'));
directiony = hdr.details.acquisition.attvalue{num};

num = find(niak_find_str_cell(hdr.details.acquisition.varatts,'direction_z'));
directionz = hdr.details.acquisition.attvalue{num};

%size(bvalues,1)
	%					      size(directionx,1)
		%				      size(directiony,1)
			%			      size(directionz,1)
% have to handle when baseline bvalue is not exactly 0
b_baseline = min(bvalues);

for i=1:size(bvalues,2)
  if (bvalues(i)==b_baseline)
    zeroimage=i;
    ab(i)=0;
  end
  if (bvalues(i)~=b_baseline)
    ab(i)=(bvalues(i)/(directionx(i)^2+directiony(i)^2+directionz(i)^2));
  end
  bxx(i)=directionx(i)*directionx(i)*ab(i);
  byy(i)=directiony(i)*directiony(i)*ab(i);
  bzz(i)=directionz(i)*directionz(i)*ab(i);
  bxy(i)=directionx(i)*directiony(i)*ab(i);
  byz(i)=directiony(i)*directionz(i)*ab(i);
  bxz(i)=directionx(i)*directionz(i)*ab(i);
  
  
end

number_images=size(bxx,2);
b=[bxx' byy' bzz' bxy' bxz' byz'];

    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hh = vectorplot(varargin)

% modified from quiver.m: s is now a matrix so scaling is different for each
% arrow.  No arrowheads (but option is left in code).


% Arrow head parameters
alpha = 0.33; 	% Size of arrow head relative to the length of the vector
beta = 0.33;  % Width of the base of the arrow head relative to the length
autoscale = 1; % Autoscale if ~= 0 then scale by this.
plotarrows = 0; % Plot arrows
sym = '';

filled = 0;
ls = '-';
ms = '';
col = 'k';

nin = nargin;
% Parse the string inputs
while ischar(varargin{nin}),
  vv = varargin{nin};
  if ~isempty(vv) & strcmp(lower(vv(1)),'f')
    filled = 1;
    nin = nin-1;
  else
    [l,c,m,msg] = colstyle(vv);
    if ~isempty(msg), 
      error(sprintf('Unknown option "%s".',vv));
    end
    if ~isempty(l), ls = l; end
    if ~isempty(c), col = c; end
    if ~isempty(m), ms = m; plotarrows = 0; end
    if isequal(m,'.'), ms = ''; end % Don't plot '.'
    nin = nin-1;
  end
end

error(nargchk(2,5,nin));

% Check numeric input arguments
if nin<4, % quiver(u,v) or quiver(u,v,s)
  [msg,x,y,u,v] = xyzchk(varargin{1:2});
else
  [msg,x,y,u,v] = xyzchk(varargin{1:4});
end
if ~isempty(msg), error(msg); end

if nin==3 | nin==5, % quiver(u,v,s) or quiver(x,y,u,v,s)  - s is matrix 
  scalematrix = varargin{nin};
end


% expand u,v if they are scalars (they aren't for vectorplot.m...)
if prod(size(u))==1, u = u(ones(size(x))); end
if prod(size(v))==1, v = v(ones(size(u))); end

if autoscale,
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n;
  dely = diff([min(y(:)) max(y(:))])/m;
  %len = sqrt((u.^2 + v.^2)/(delx.^2 + dely.^2));
  as=ones(size(u));
  for i=1:size(u,1)
    for j=1:size(u,2)
      as(i,j) = (scalematrix(i,j)); 
      u(i,j) = u(i,j)*as(i,j); v(i,j) = v(i,j)*as(i,j);
    end;
  end;
end


ax = newplot;
next = lower(get(ax,'NextPlot'));
hold_state = ishold;

% Make velocity vectors
x = x(:).'; y = y(:).';
u = u(:).'; v = v(:).';
uu = [x;x+u;repmat(NaN,size(u))];
vv = [y;y+v;repmat(NaN,size(u))];

h1 = plot(uu(:),vv(:),[col ls]); %uu(:) regards uu as a column

if plotarrows,
  % Make arrow heads and plot them
  hu = [x+u-alpha*(u+beta*(v+eps));x+u; ...
        x+u-alpha*(u-beta*(v+eps));repmat(NaN,size(u))];
  hv = [y+v-alpha*(v-beta*(u+eps));y+v; ...
        y+v-alpha*(v+beta*(u+eps));repmat(NaN,size(v))];
  hold on
  h2 = plot(hu(:),hv(:),[col ls]);
else
  h2 = [];
end

if ~isempty(ms), % Plot marker on base
  hu = x; hv = y;
  hold on
  h3 = plot(hu(:),hv(:),[col ms]);
  if filled, set(h3,'markerfacecolor',get(h1,'color')); end
else
  h3 = [];
end

if ~hold_state, hold off, view(2); set(ax,'NextPlot',next); end

if nargout>0, hh = [h1;h2;h3]; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

