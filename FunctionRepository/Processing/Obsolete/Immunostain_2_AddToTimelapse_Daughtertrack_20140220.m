function Immunostain_2_AddToTimelapse(row,col,site)
%row=4;col=8;site=1;
row=5;col=2;site=1;
%row='E';col='04';site='1';
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
%experimentpath='Kyuho\';
experimentpath='2014-02-08_H2B-DHB-p21dCy1dK\';

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
name1='H2B_'; %nuc
nucedgename='nucedge_';
name2='DHB_';
name3='p21dCy1dK_';
IFname1='p21_';
IFname2='DAPI_';
moviebin=1;
if moviebin==1
    nucr=12;
    debrisarea=100; %MCF-10A: 200
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea;
blobthreshold=-0.02;
timetotal=tic;
%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawprev=single(imread([rawdir,name1,num2str(totalframes),'.tif']));
[nuc_mask_prev,~]=blobdetector_foreground(log(rawprev),nucr,blobthreshold,debrisarea);
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'stain.tif']));
raw2=single(imread([rawdir,name2,'stain.tif']));
raw3=single(imread([rawdir,name3,'stain.tif']));
IFraw1=single(imread([rawdir,IFname1,'stain.tif']));
%IFraw2=single(imread([rawdir,IFname2,'poststain.tif']));
NEfile=[maskdir,nucedgename,'stain.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
[nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
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
fg3=raw3; fg3(nanmask)=NaN;
bg1=blocksmooth(fg1,10);
bg2=blocksmooth_mode(fg2,1);
bg3=blocksmooth(fg3,10);
real1=raw1-bg1;
real2=raw2-bg2;
real3=raw3-bg3;
%%% calculate IF background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFblur1=imfilter(IFraw1,fspecial('disk',nucr),'symmetric');
iqrmult=4;
IFhighfliers1=maskhighfliers_nonnuc(IFblur1,nanmask,nucr,iqrmult);
IFfg1=nanmask | IFhighfliers1;
IFfg=IFfg1;
%IFfg=IFfg1 | IFfg2;
IFblur1(IFfg)=NaN;
IFbg1=blocksmooth(IFblur1,10);
IFreal1=IFraw1-IFbg1;
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
%[tracked,nuc_label]=adaptivetrack_IF(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr);
debugpackage={rawprev,nuc_mask_prev,nuc_mask,prevjitter,reljitter};
[tracked,nuc_label]=adaptivetrack_IF_daughtertrack_secondtry(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,debugpackage);
numcells=size(tracked,1);
nuc_center=tracked(:,[1 2]);
nuc_area=tracked(:,3);
nuc_mass=tracked(:,4);
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_label,'remove');
if maskwrite
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellid=find(~isnan(nuc_mass));
nuc_info=regionprops(nuc_label,'PixelIdxList');
ring_label=getcytoring(nuc_label);
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN; sig2=nanvec; sig3=nanvec; sig2ring=nanvec;
IFsig1=nanvec; IFsig1ring=nanvec; IFsig2=nanvec;
for cc=cellid'
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig2ringall=real2(ring_info(cc).PixelIdxList);
     sig2ring(cc)=mean(sig2ringall(sig2ringall>=median(sig2ringall)));
    sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    IFsig1(cc)=mean(IFraw1(nuc_info(cc).PixelIdxList));
    %IFsig1(cc)=max(IFreal1(nuc_info(cc).PixelIdxList)); %use max for EdU
    %IFsig1ringall=IFreal1(ring_info(cc).PixelIdxList);
     %IFsig1ring(cc)=mean(IFsig1ringall(IFsig1ringall>=median(IFsig1ringall)));
    %IFsig2(cc)=median(IFreal2(nuc_info(cc).PixelIdxList));
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
IF_bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
IFsig1=IFsig1-IF_bleedthrough;
%IF_sig_max=IF_sig_max-IF_bleedthrough;
%IF_sig_mean=IF_sig_mean-IF_bleedthrough;
%%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
%load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
%nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
%DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
%IFsig2=IFsig2-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,IFsig1,IFsig1ring,IFsig2];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,IFsig1];
save([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(IFraw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}