function Immunostain_2_AddToTimelapse(row,col,site)
row=2;col=2;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
experimentpath='Kyuho\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data\']);
%rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
%maskdir = [imagepath,experimentpath,'Mask\',shot];
maskdir = [imagepath,experimentpath,'Mask'];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='CFP_'; %nuc
nucedgename='nucedge_';
name2='TexasRed_';
name3='YFP_';
IFname1='pS6_';
IFname2='Hoechst_';
moviebin=1;
if moviebin==1
    nucr=12;
    debrisarea=200; %MCF-10A: 200
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
blobthreshold=-0.02;
timetotal=tic;
%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask_prev=single(imread([maskdir,'\',nucedgename,num2str(totalframes),'.tif']));
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'stain.tif']));
raw2=single(imread([rawdir,name2,'stain.tif']));
raw3=single(imread([rawdir,name3,'stain.tif']));
IFraw1=single(imread([rawdir,IFname1,'stain.tif']));
IFraw2=single(imread([rawdir,IFname2,'stain.tif']));
NEfile=[maskdir,'\',nucedgename,'stain.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask_withborders=nuc_mask;
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label,numcells]=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,raw1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_bg=centroidbackground(nuc_mask_withborders,nuc_center,raw1,nucr,3,0.25,50);
nuc_density=nuc_density-nuc_bg;
%%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mass=nuc_density.*nuc_area;
%%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx,reljity]=detectGlobalShift(nuc_mask_prev,nuc_mask);
reljitx=-reljitx;
reljity=-reljity;
IFjitter=jitters(totalframes,:)+[reljitx,reljity];
nuc_center(:,1)=nuc_center(:,1)+IFjitter(1);
nuc_center(:,2)=nuc_center(:,2)+IFjitter(2);
%%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
%%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
areaonly=1;
[tracked,nuc_label]=adaptivetrack_IF(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,areaonly);
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
nanvec=ones(numcells,1)*NaN; sig2=nanvec; sig3=nanvec; sig2ring=nanvec;
IFsig1=nanvec; IFsig1ring=nanvec; IFsig2=nanvec;
for cc=cellid'
    sig2(cc)=median(raw2(nuc_info(cc).PixelIdxList));
    sig2ringall=raw2(ring_info(cc).PixelIdxList);
     sig2ring(cc)=mean(sig2ringall(sig2ringall>=median(sig2ringall)));
    sig3(cc)=median(raw3(nuc_info(cc).PixelIdxList));
    IFsig1(cc)=median(IFraw1(nuc_info(cc).PixelIdxList));
    IFsig1ringall=IFraw1(ring_info(cc).PixelIdxList);
     IFsig1ring(cc)=mean(IFsig1ringall(IFsig1ringall>=median(IFsig1ringall)));
    IFsig2(cc)=median(IFraw2(nuc_info(cc).PixelIdxList));
end
%%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bg2=centroidbackground(nuc_mask_withborders,nuc_center(cellid,:),raw2,nucr,10,0.25,10);
bg3=centroidbackground(nuc_mask_withborders,nuc_center(cellid,:),raw3,nucr,3,0.25,50);
IFbg1=centroidbackground(nuc_mask_withborders,nuc_center(cellid,:),IFraw1,nucr,10,0.25,10);
IFbg2=centroidbackground(nuc_mask_withborders,nuc_center(cellid,:),IFraw2,nucr,3,0.25,50);
sig2(cellid)=sig2(cellid)-bg2;
sig2ring(cellid)=sig2ring(cellid)-bg2;
sig3(cellid)=sig3(cellid)-bg3;
IFsig1(cellid)=IFsig1(cellid)-IFbg1;
IFsig1ring(cellid)=IFsig1ring(cellid)-IFbg1;
IFsig2(cellid)=IFsig2(cellid)-IFbg2;
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
% IF_bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
% IFsig1=IFsig1-IF_bleedthrough;
% IF_sig_max=IF_sig_max-IF_bleedthrough;
% IF_sig_mean=IF_sig_mean-IF_bleedthrough;
% %%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
% nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
% DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
% IFsig2=IFsig2-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,IFsig1,IFsig1ring,IFsig2];
save([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(IFraw2));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}