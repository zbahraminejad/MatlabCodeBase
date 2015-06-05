function DefineObjects
%%% DESCRIPTION
% This is a short script to identify foreground objects. I have already
% subtracted the background fluorescence, so this looks better than the raw
% image coming right off of the microscope. The identified objects will
% often consist of multiple nuclei, which illustrates the insufficiency of
% the code as is.

%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawdir = 'H:\Documents\Projects\Keshav\FixedImage\';

%%% constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuclearradius=12; %avg nuclear radius
debrisarea=200; %min# pixels to be considered a cell
blobthreshold=-0.03; %parameter

%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawmask=single(imread([rawdir,'Hoechst.tif']));
rawsensor=single(imread([rawdir,'DHB.tif']));

%%% detect contiguous objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma=0.75*nuclearradius/sqrt(2); %set parameter for filter kernel
h=sigma^2*fspecial('log',[nuclearradius*2 nuclearradius*2],sigma); %define filter kernel
nuc_log=imfilter(rawmask,h,'symmetric'); %filter image
nuc_mask=nuc_log<blobthreshold; %define mask by thresholding on filtered image
nuc_mask=imfill(nuc_mask,'holes'); %remove internal holes
nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %remove connections
nuc_mask=~bwmorph(~nuc_mask,'diag'); %remove connections
nuc_mask=~bwmorph(~nuc_mask,'bridge'); %remove connections
nuc_mask=bwareaopen(nuc_mask,debrisarea); %remove objects that fall below min# pixel threshold

%%% view images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(rawmask));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);