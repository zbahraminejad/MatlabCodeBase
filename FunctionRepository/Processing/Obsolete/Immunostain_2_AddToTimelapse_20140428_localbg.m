function Immunostain_2_AddToTimelapse(row,col,site)
%row=1;col=8;site=1;
%row='E';col='04';site='1';
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
shadingpath='H:\Images\ShadingImages\20140410 DCYTC 20x\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
%experimentpath='Kyuho\';
%experimentpath='2014-02-08_H2B-DHB-p21dCy1dK\';
%experimentpath='2014-03-18_CyclingCDK4iWee1iCDC25i\';
%experimentpath='20131213 R-point CC\20140405 2C siRNA CDK4i\';
experimentpath='20131213 R-point CC\20140407 2C CDK2 Hysteresis\';
%experimentpath='20131213 R-point CC\20140413 2C DHB Hysteresis\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=([projectpath,experimentpath,'Data\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask',shot,'_'];
end
maskwrite=0;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='H2B_';
nucedgename='nucedge_';
name2='DHB_';
%name3='pRb_';
IFname1='tRb_';
IFname2='pRb_';
%IFname3='EdU_';
moviebin=1;
if moviebin==1
    %nucr=12;
    %debrisarea=100; %MCF-10A: 100
    nucr=20;
    debrisarea=600; %MCF-10A 20xBin1
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea;
%blobthreshold=-0.02; %MCF10A 10xBin1
blobthreshold=-0.03; %MCF10A 20xBin1
timetotal=tic;
%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawprev=single(imread([rawdir,name1,num2str(totalframes),'.tif']));
[nuc_mask_prev,~]=blobdetector_foreground(log(rawprev),nucr,blobthreshold,debrisarea);
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
% load([shadingpath,'CFP','.mat'],'shadingcorrection'); sc1=shadingcorrection;
% load([shadingpath,'YFP','.mat'],'shadingcorrection'); sc2=shadingcorrection;
% load([shadingpath,'DAPI','.mat'],'shadingcorrection'); sc3=shadingcorrection;
% load([shadingpath,'TxRed','.mat'],'shadingcorrection'); sc4=shadingcorrection;
% load([shadingpath,'Cy5','.mat'],'shadingcorrection'); sc5=shadingcorrection;
load([shadingpath,'BG_',num2str(site),'.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([shadingpath,'CFP_',num2str(site),'.mat'],'shadingcorrection'); sc1=shadingcorrection;
load([shadingpath,'YFP_',num2str(site),'.mat'],'shadingcorrection'); sc2=shadingcorrection;
load([shadingpath,'TxRed_',num2str(site),'.mat'],'shadingcorrection'); sc3=shadingcorrection;
load([shadingpath,'Cy5_',num2str(site),'.mat'],'shadingcorrection'); sc4=shadingcorrection;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'stain.tif'])); raw1=(raw1-bgcmos)./sc1;
raw2=single(imread([rawdir,name2,'stain.tif'])); raw2=(raw2-bgcmos)./sc2;
%raw3=single(imread([rawdir,name3,'stain.tif']));
IFraw1=single(imread([rawdir,IFname1,'stain.tif'])); IFraw1=(IFraw1-bgcmos)./sc3;
IFraw2=single(imread([rawdir,IFname2,'stain.tif'])); IFraw2=(IFraw2-bgcmos)./sc4;
%IFraw3=single(imread([rawdir,IFname3,'stain.tif'])); IFraw3=(IFraw3-bgcmos)./sc5;
NEfile=[maskdir,nucedgename,'stain.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
[nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% sensor background subtract %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanmask=imdilate(foreground,strel('disk',nucr));
nanmaskcyto=imdilate(foreground,strel('disk',2*nucr));
fg1=raw1; fg1(nanmask)=NaN;
fg2=raw2; fg2(nanmaskcyto)=NaN;
%fg3=raw3; fg3(nanmask)=NaN;
bg1=blocksmooth(fg1,10);
bg2=blocksmooth_mode(fg2,1);
%bg3=blocksmooth(fg3,10);
real1=raw1-bg1;
real2=raw2-bg2;
%real3=raw3-bg3;
%%% calculate IF background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFblur1=imfilter(IFraw1,fspecial('disk',nucr),'symmetric');
IFblur2=imfilter(IFraw2,fspecial('disk',nucr),'symmetric');
%IFblur3=imfilter(IFraw3,fspecial('disk',nucr),'symmetric');
iqrmult=4;
IFhighfliers1=maskhighfliers_nonnuc(IFblur1,nanmask,nucr,iqrmult);
IFhighfliers2=maskhighfliers_nonnuc(IFblur2,nanmask,nucr,iqrmult);
%IFhighfliers3=maskhighfliers_nonnuc(IFblur3,nanmask,nucr,iqrmult);
IFfg1=nanmask | IFhighfliers1;
IFfg2=nanmask | IFhighfliers2;
%IFfg3=nanmask | IFhighfliers3;
%IFfg=IFfg1;
IFfg=IFfg1 | IFfg2;
%IFfg=IFfg1 | IFfg2 | IFfg3;
IFblur1(IFfg)=NaN;
IFblur2(IFfg)=NaN;
%IFblur3(IFfg)=NaN;
IFbg1=blocksmooth(IFblur1,10);
IFbg2=blocksmooth(IFblur2,10);
%IFbg3=blocksmooth(IFblur3,10);
IFreal1=IFraw1-IFbg1;
IFreal2=IFraw2-IFbg2;
%IFreal3=IFraw3-IFbg3;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label,numcells]=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mass=nuc_density.*nuc_area;
%%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx,reljity]=detectGlobalShift(nuc_mask_prev,nuc_mask);
reljitx=-reljitx;
reljity=-reljity;
reljitter=[reljitx,reljity];
prevjitter=jitters(totalframes,:);
IFjitter=prevjitter+reljitter;
nuc_center(:,1)=nuc_center(:,1)+IFjitter(1);
nuc_center(:,2)=nuc_center(:,2)+IFjitter(2);
%%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
%%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
debugpackage={rawprev,nuc_mask_prev,nuc_mask,prevjitter,reljitter};
[tracked,nuc_label]=adaptivetrack_IF(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,debugpackage);
numcells=size(tracked,1);
nuc_center=tracked(:,[1 2]);
nuc_area=tracked(:,3);
nuc_mass=tracked(:,4);
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if maskwrite
    extractmask=bwmorph(nuc_label,'remove');
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellid=find(~isnan(nuc_mass));
nuc_info=regionprops(nuc_label,'PixelIdxList');
ring_label=getcytoring(nuc_label);
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec; sig2ring=nanvec;
IFsig1=nanvec; IFsig2=nanvec; IFsig3=nanvec;
for cc=cellid'
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig2ringall=real2(ring_info(cc).PixelIdxList);
     sig2ring(cc)=median(sig2ringall(sig2ringall>=median(sig2ringall)));
    %sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    IFsig1(cc)=mean(IFreal1(nuc_info(cc).PixelIdxList));
    IFsig2(cc)=mean(IFreal2(nuc_info(cc).PixelIdxList));
    %IFsig3(cc)=mean(IFreal3(nuc_info(cc).PixelIdxList));
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
%IF_bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
%IFsig1=IFsig1-IF_bleedthrough;
%IF_sig_max=IF_sig_max-IF_bleedthrough;
%IF_sig_mean=IF_sig_mean-IF_bleedthrough;
%%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
%load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
%nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
%DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
%IFsig2=IFsig2-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,IFsig1,IFsig1ring,IFsig2];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig2ring,IFsig1];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring,IFsig1,IFsig2];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring,IFsig1,IFsig2,IFsig3];
save([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}