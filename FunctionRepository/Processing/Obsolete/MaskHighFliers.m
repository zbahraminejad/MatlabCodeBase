function Immunostain_2(row,col,site)
row=2;col=1;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131217\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='DAPI_'; %nuc
name2='FITC_';
name3='Cy3_';
name4='Cy5_';
timetotal=tic;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%raw1=single(imread([rawdir,name1,'stain.tif']));
%raw2=single(imread([rawdir,name2,'stain.tif']));
raw3=single(imread([rawdir,name3,'stain.tif']));
%raw4=single(imread([rawdir,name4,'stain.tif']));
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=150;
debrisarea=round(pi*100^2);
blobthreshold=-0.02;
%nuc_mask=smeardetector(log(raw3),nucr,blobthreshold,debrisarea);
nuc_mask=smeardetector(log(imfilter(raw3,fspecial('disk',150),'symmetric')),nucr,blobthreshold,debrisarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask_withborder=nuc_mask;
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw3));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}