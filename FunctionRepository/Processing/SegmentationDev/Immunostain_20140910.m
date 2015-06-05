%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawdir='D:\Documents\MATLAB\Imaging\Processing\SegmentationDev\';
%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='CaseMCF10A_Hoechst';
nucr=12;
debrisarea=200;
%nucr=20;
%debrisarea=700;
boulderarea=20*debrisarea;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'.tif']));

%%%%%%%%%%%%%%%% ALTERNATIVE SEGMENTATION APPROACHES %%%%%%%%%%%%%%%%%%%%%%

%%% Segmentdeflections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=threshmask(raw1,3);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
%nuc_mask=segmentdeflections_leavedebris(nuc_mask,nucr);
nuc_mask=segmentdeflections_bwboundaries_leavedebris(nuc_mask,nucr);
keyboard;

%{
%%% Laplacian of Gaussian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blobthreshold=-0.03;
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground=nuc_mask;
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%}

%{
%%% marker-based watershed: erosion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=threshmask(raw1,3);
nuc_mask=markershed(nuc_mask,10);
%nuc_mask=markershed_regionalmax(nuc_mask,raw1);
foreground=nuc_mask;
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_leavedebris(nuc_mask,nucr);
%nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
%nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%hist(nuc_area,100);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
%}

%{
%%% marker-based watershed: regionalmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=threshmask(raw1);
nuc_mask=markershed_regionalmax(nuc_mask,raw1,nucr);
foreground=nuc_mask;
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_leavedebris(nuc_mask,nucr);
%nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
%nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%hist(nuc_area,100);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
%}

%{
%%% screen by perimeter/area *Didn't screen anything out %%%%%%%%%%%%%%%%%%
nuc_info=struct2cell(regionprops(nuc_mask,'Area','Perimeter')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_perim=squeeze(cell2mat(nuc_info(2,1,:)));
%hist(nuc_perim,100); debrisperim=90;
perimarearatio=nuc_perim./nuc_area;
%hist(perimarearatio,100);
perimareathresh=0.13;
unsegidx=find(perimarearatio>perimareathresh);
nuc_pix=regionprops(nuc_mask,'PixelIdxList');
for i=unsegidx'
    nuc_mask(nuc_pix(i).PixelIdxList)=0;
end
%}

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}