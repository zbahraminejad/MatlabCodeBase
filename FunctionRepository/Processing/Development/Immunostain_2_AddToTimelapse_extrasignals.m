function Immunostain_2_AddToTimelapse(row,col,site)
%row='B';col='07';site='1';
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data\']);
rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
maskdir = [imagepath,experimentpath,'Mask\',shot];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12; %calculate avg nuclear radius with Measure_avgnucrad.m
nucname='CFP_'; %nuc
nucedgename='nucedge_';
YFPname='YFP_';
RFPname='TexasRed_';
IFname='Cy5_';
DAPIname='DAPI_';
blobthreshold=-0.02;
timetotal=tic;
%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_raw_prev=single(imread([rawdir,nucname,num2str(totalframes),'.tif']));
nuc_mask_prev=blobdetector(log(nuc_raw_prev),nucr,blobthreshold);
[height,width]=size(nuc_mask_prev);
nuc_mask_prev([1 height],1:width)=1; nuc_mask_prev(1:height,[1 width])=1;
border=bwareaopen(nuc_mask_prev,height*2+width*2);
nuc_mask_prev=logical(nuc_mask_prev-border);
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_raw=single(imread([rawdir,nucname,'poststain.tif']));
YFP_raw=single(imread([rawdir,YFPname,'poststain.tif']));
RFP_raw=single(imread([rawdir,RFPname,'poststain.tif']));
IF_raw=single(imread([rawdir,IFname,'poststain.tif']));
DAPI_raw=single(imread([rawdir,DAPIname,'poststain.tif']));
NEfile=[maskdir,'\',nucedgename,'poststain.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=blobdetector(log(nuc_raw),nucr,blobthreshold);
%nuc_mask=segmentdeflections(nuc_mask,nucr,0); %attempt on all cells
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2);
nuc_mask=logical(nuc_mask-border);
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,nuc_raw,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_bg=getbackground_H2B(nuc_mask,nuc_center,nuc_raw,nucr,0.25);
nuc_density=nuc_density-nuc_bg;
nuc_mass=nuc_density.*nuc_area;
%%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx,reljity]=detectGlobalShift(imfill(nuc_mask_prev,'holes'),nuc_mask);
reljitx=-reljitx;
reljity=-reljity;
IFjitter=jitters(totalframes,:)+[reljitx,reljity];
nuc_center(:,1)=nuc_center(:,1)+IFjitter(1);
nuc_center(:,2)=nuc_center(:,2)+IFjitter(2);
%%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
%%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
[tracked,nuc_label]=adaptivetrack_IF_area(tracedata(:,totalframes,1:3),initdata,nuc_label,nucr);
numcells=size(tracked,1);
nuc_center=tracked(:,[1 2]);
nuc_area=tracked(:,3);
nuc_mass=tracked(:,4);
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_label,'remove');
imwrite(uint16(extractmask),NEfile);
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellid=find(~isnan(nuc_mass));
nuc_info=regionprops(nuc_label,'PixelIdxList');
ring_label=getcytoring(nuc_label);
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN; YFP_sig=nanvec; RFP_sig=nanvec; YFP_sigring=nanvec;
DAPI_sig=nanvec; IF_sig_med=nanvec; sig_max=nanvec; sig_mean=nanvec;
for cc=cellid'
    YFP_sig(cc)=median(YFP_raw(nuc_info(cc).PixelIdxList));
    RFP_sig(cc)=median(RFP_raw(nuc_info(cc).PixelIdxList));
    ringall=YFP_raw(ring_info(cc).PixelIdxList);
     YFP_sigring(cc)=mean(ringall(ringall>=median(ringall)));
    DAPI_sig(cc)=mean(DAPI_raw(nuc_info(cc).PixelIdxList));
    IF_sig_med(cc)=median(IF_raw(nuc_info(cc).PixelIdxList));
    sig_max(cc)=max(RFP_raw(nuc_info(cc).PixelIdxList));
    sig_mean(cc)=mean(RFP_raw(nuc_info(cc).PixelIdxList));
end
%%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[YFP_bg,RFP_bg]=getbackground_12(nuc_mask,nuc_center,YFP_raw,RFP_raw,nucr,20,0.25,10,50);
YFP_sig(cellid)=YFP_sig(cellid)-YFP_bg;
RFP_sig(cellid)=RFP_sig(cellid)-RFP_bg;
YFP_sigring(cellid)=YFP_sigring(cellid)-YFP_bg;
[DAPI_bg,IF_bg]=getbackground_12(nuc_mask,nuc_center,DAPI_raw,IF_raw,nucr,10,0.25,50,50);
DAPI_sig(cellid)=DAPI_sig(cellid)-DAPI_bg;
IF_sig_med(cellid)=IF_sig_med(cellid)-IF_bg;
sig_max(cellid)=sig_max(cellid)-RFP_bg;
sig_mean(cellid)=sig_mean(cellid)-RFP_bg;
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
IF_bleedthrough=bleedthroughrate(1)+RFP_sig*bleedthroughrate(2);
IF_sig_med=IF_sig_med-IF_bleedthrough;
%sig_max=sig_max-IF_bleedthrough;
%sig_mean=sig_mean-IF_bleedthrough;
%%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
DAPI_sig=DAPI_sig-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,YFP_sig,RFP_sig,YFP_sigring,DAPI_sig,IF_sig_med,sig_max,sig_mean];
save([datadir,'sensorIF_',shot,'.mat'],'IFdata','IFjitter');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}