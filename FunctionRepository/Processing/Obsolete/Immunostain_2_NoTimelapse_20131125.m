function Immunostain_2_NoTimelapse(row,col,site)
row=1;col=2;site=2;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath='2013-11-16_CDK4_pRb_Abs\Experiment_2013-11-16_S780\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data\']);
rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
maskdir = [imagepath,experimentpath,'Mask\',shot];
if ~exist(maskdir,'dir')
    mkdir(maskdir);
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucname='Hoechst_'; %nuc
YFPname='FITC_';
RFPname='Cy5_';
nucedgename='nucedge_';
moviebin=1;
if moviebin==1
    nucr=12;
    debrisarea=200; %MCF-10A: 200
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
blobthreshold=-0.03;
timetotal=tic;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_raw=single(imread([rawdir,nucname,'1.tif']));
YFP_raw=single(imread([rawdir,YFPname,'1.tif']));
RFP_raw=single(imread([rawdir,RFPname,'1.tif']));
NEfile=[maskdir,'\',nucedgename,'1.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%YFP_mask=smallblobdetector(log(YFP_raw),12,-0.03,20,200);
nuc_mask=blobdetector(log(nuc_raw),nucr,blobthreshold,debrisarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_raw);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=struct2cell(regionprops(nuc_mask,nuc_raw,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_bg=getbackground_H2B(nuc_mask,nuc_center,nuc_raw,nucr,0.25);
nuc_density=nuc_density-nuc_bg;
nuc_mass=nuc_density.*nuc_area;
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
imwrite(uint16(extractmask),NEfile);
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=numel(nuc_mass);
nuc_info=regionprops(nuc_mask,'PixelIdxList');
%ring_label=getcytoring(nuc_label);
%ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN; YFP_sig=nanvec; RFP_sig=nanvec; YFP_sigring=nanvec;
YFP_sig_mean=nanvec; YFP_sig_max=nanvec; RFP_sig_mean=nanvec; RFP_sig_max=nanvec;
for cc=1:numcells
    YFP_sig(cc)=median(YFP_raw(nuc_info(cc).PixelIdxList));
    YFP_sig_mean(cc)=mean(YFP_raw(nuc_info(cc).PixelIdxList));
    YFP_sig_max(cc)=max(YFP_raw(nuc_info(cc).PixelIdxList));
    RFP_sig(cc)=median(RFP_raw(nuc_info(cc).PixelIdxList));
    RFP_sig_mean(cc)=mean(RFP_raw(nuc_info(cc).PixelIdxList));
    RFP_sig_max(cc)=max(RFP_raw(nuc_info(cc).PixelIdxList));
    %ringall=YFP_raw(ring_info(cc).PixelIdxList);
     %YFP_sigring(cc)=mean(ringall(ringall>=median(ringall)));
end
%%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[YFP_bg,RFP_bg]=getbackground_12(nuc_mask,nuc_center,YFP_raw,RFP_raw,nucr,20,0.25,10,50);
YFP_sig=YFP_sig-YFP_bg;
YFP_sig_mean=YFP_sig_mean-YFP_bg;
YFP_sig_max=YFP_sig_max-YFP_bg;
RFP_sig=RFP_sig-RFP_bg;
RFP_sig_mean=RFP_sig_mean-RFP_bg;
RFP_sig_max=RFP_sig_max-RFP_bg;
%YFP_sigring(cellid)=YFP_sigring(cellid)-YFP_bg;
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
% IF_bleedthrough=bleedthroughrate(1)+RFP_sig*bleedthroughrate(2);
% IF_sig_med=IF_sig_med-IF_bleedthrough;
% IF_sig_max=IF_sig_max-IF_bleedthrough;
% IF_sig_mean=IF_sig_mean-IF_bleedthrough;
%%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
% nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
% DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
% DAPI_sig=DAPI_sig-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,YFP_sig,RFP_sig,YFP_sig_mean,RFP_sig_mean,YFP_sig_max,RFP_sig_max];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempimage=mat2gray(YFP_raw);
tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}