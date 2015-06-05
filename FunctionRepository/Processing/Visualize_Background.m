function Timelapse(row,col,site)
row=1;col=4;site=1;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
%shadingpath= 'H:\Images\ShadingImages\20140402_DAPI_CFP_YFP_TxRed_Cy5\';
shadingpath='H:\Images\ShadingImages\20140410 DCYTC 20x\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
experimentpath='20131213 R-point CC\20140407 2C CDK2 Hysteresis\';
%experimentpath='Heewon\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=([projectpath,experimentpath,'Data\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,shot,'\',shot,'_']; %Heewon
    maskdir=[imagepath,experimentpath,'Mask\',shot];
    %maskdir=[imagepath,experimentpath,'Mask\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask'];
end
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
if separatedirectories==0
    maskdir=[maskdir,'\',shot,'_'];
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
moviebin=1;
if moviebin==1
    %nucr=12;
    %debrisarea=100; %MCF-10A 10xBin1: H2B:100 NLS:200
    nucr=20;
    debrisarea=600;
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
%blobthreshold=-0.02;  %default 10xbin1=-0.02
blobthreshold=-0.03;
name1='H2B_'; %nuc
name2='DHB_'; %DHB
f=1;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load([shadingpath,'CFP','.mat'],'shadingcorrection'); sc1=shadingcorrection;
%load([shadingpath,'YFP','.mat'],'shadingcorrection'); sc2=shadingcorrection;

load([shadingpath,'BG_',num2str(site),'.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([shadingpath,'CFP_',num2str(site),'.mat'],'shadingcorrection'); sc1=shadingcorrection;
load([shadingpath,'YFP_',num2str(site),'.mat'],'shadingcorrection'); sc2=shadingcorrection;

raw1=single(imread([rawdir,name1,num2str(f),'.tif']));
raw1=raw1-bgcmos;
raw1=raw1./sc1;
raw2=single(imread([rawdir,name2,num2str(f),'.tif'])); %raw2=raw2./sc2;
[~,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nanmask=imdilate(foreground,strel('disk',nucr));
nanmaskcyto=imdilate(foreground,strel('disk',2*nucr));
fg1=raw1; fg1(nanmask)=NaN;
fg2=raw2; fg2(nanmaskcyto)=NaN;
bg1=blocksmooth(fg1,10);
bg2=blocksmooth(fg2,10);
keyboard;
imagesc(bg2);colorbar;
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=imadjust(mat2gray(raw1));
tempframe(:,:,3)=0;
tempframe(:,:,1)=extractmask;
figure,imshow(tempframe);
%}