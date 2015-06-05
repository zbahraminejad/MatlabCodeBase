function Immunostain_2_MSC_bgsubfirst_sparse(row,col,site)
%row=2;col=4;site=16;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
%imagepath = 'E:\';
imagepath = 'H:\Images\';
%shadingpath='H:\Images\ShadingImages\20140411_DFCTC_10x\';
shadingpath='H:\Images\ShadingImages\20140410 DCYTC 20x\';
%experimentpath='20131213 R-point CC\20140406 siRNA CDK4i\';
experimentpath='20131126 R-point SR\20140415 CDK46hysteresis-under-CDK2i\';
%experimentpath='20131213 R-point CC\20140413 2C DHB Hysteresis\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%datadir = ([projectpath,experimentpath,'Data\']);
datadir = ([projectpath,experimentpath,'Data_min\']);
separatedirectories=0;
if separatedirectories
    rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir = [imagepath,experimentpath,'Mask\',shot];
else
    rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir = [imagepath,experimentpath,'Mask'];
end
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='Hoechst'; %nuc
name2='DHB';
name3='tRb';
name4='pRb';
nucedgename='nucedge';
moviebin=1;
if moviebin==1
    nucr=12; %MCF-10A:12 YT+:8
    debrisarea=200; %MCF-10A:200 YT+:300
    %nucr=20;
    %debrisarea=600;
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=10*debrisarea; %MCF-10A:2000 YT+:2000
blobthreshold=-0.03; %MCF-10A:-0.03 YT+:-0.03
timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG_bin2','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([shadingpath,'DAPI_bin2','.mat'],'shadingcorrection'); sc1=shadingcorrection;
%load([shadingpath,'CFP','.mat'],'shadingcorrection'); sc1=shadingcorrection;
%load([shadingpath,'FITC','.mat'],'shadingcorrection'); sc2=shadingcorrection;
load([shadingpath,'YFP_bin2','.mat'],'shadingcorrection'); sc2=shadingcorrection;
load([shadingpath,'TxRed_bin2','.mat'],'shadingcorrection'); sc3=shadingcorrection;
load([shadingpath,'Cy5_bin2','.mat'],'shadingcorrection'); sc4=shadingcorrection;

% load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
% load([shadingpath,'DAPI','.mat'],'shadingcorrection'); sc1=shadingcorrection;
% load([shadingpath,'YFP','.mat'],'shadingcorrection'); sc2=shadingcorrection;
% load([shadingpath,'TxRed','.mat'],'shadingcorrection'); sc3=shadingcorrection;
% load([shadingpath,'Cy5','.mat'],'shadingcorrection'); sc4=shadingcorrection;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif'])); raw1=(raw1-bgcmos)./sc1;
raw2=single(imread([rawdir,name2,'_stain.tif'])); raw2=(raw2-bgcmos)./sc2;
raw3=single(imread([rawdir,name3,'_stain.tif'])); raw3=(raw3-bgcmos)./sc3;
raw4=single(imread([rawdir,name4,'_stain.tif'])); raw4=(raw4-bgcmos)./sc4;
if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
[nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %comment out for sparse YT analysis
%boulderarea=2000; nuc_mask=blobdetector_excludelargeandwarped(raw1,nucr,blobthreshold,debrisarea,boulderarea);
%%% remove objects that are highly eccentric %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% calculate background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blur1=imfilter(raw1,fspecial('disk',nucr),'symmetric');
blur2=imfilter(raw2,fspecial('disk',nucr),'symmetric');
blur3=imfilter(raw3,fspecial('disk',nucr),'symmetric');
blur4=imfilter(raw4,fspecial('disk',nucr),'symmetric');
c1=min(blur1(:));
c2=min(blur2(:));
c3=min(blur3(:));
c4=min(blur4(:));
bg1=ones(size(nuc_mask))*c1;
bg2=ones(size(nuc_mask))*c2;
bg3=ones(size(nuc_mask))*c3;
bg4=ones(size(nuc_mask))*c4;
real1=raw1-bg1;
real2=raw2-bg2;
real3=raw3-bg3;
real4=raw4-bg4;
%%% remove nucleoli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% feature_mask=smallblobdetector(log(raw2),12,-0.03,20,150);
% nuc_mask=logical(nuc_mask.*~feature_mask);
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label=bwlabel(nuc_mask);
%smearmask=highfliers1;
%smearmask=highfliers1 | highfliers2 | highfliers3;
%smearmask=highfliers1 | highfliers2 | highfliers3 | highfliers4;
smearmask=zeros(size(nuc_label));
smear_label=nuc_label.*smearmask;
nuc_smeared=unique(smear_label(:));
nuc_smeared(1)=[]; %remove 0
nuc_label(ismember(nuc_label,nuc_smeared))=0;
nuc_mask=nuc_label>0;
nuc_info=struct2cell(regionprops(nuc_mask,'Area','Centroid')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
if maskwrite
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=numel(nuc_area);
nuc_info=regionprops(nuc_mask,'PixelIdxList');
%ring_label=getcytoring(nuc_label);
ring_label=getcytoring_2(nuc_label,2,raw2); %20xBin1:4 20xBin2:2
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig1=nanvec;
sig2=nanvec;
sig3=nanvec;
sig4=nanvec; sig2ring_75th=nanvec; sig2ring_50th=nanvec; sig2ring_nobg=nanvec;
for cc=1:numcells
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=mean(real3(nuc_info(cc).PixelIdxList));
    sig4(cc)=mean(real4(nuc_info(cc).PixelIdxList));
    ringall=real2(ring_info(cc).PixelIdxList); %hist(ringall,20);
     sig2ring_75th(cc)=prctile(ringall,75); %previously mean
     sig2ring_50th(cc)=median(ringall);
     sig2ring_nobg(cc)=median(ringall(ringall>50)); %20xBin1:75 20xBin2:50
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
%bleedthrough=bleedthroughrate(1)+sig2*bleedthroughrate(2);
%sig3=sig3-bleedthrough;
%%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
% nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
% DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
% DAPI_sig=DAPI_sig-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig2ring_75th,sig2ring_50th,sig2ring_nobg,sig3,sig4];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4,sig2_ring,sig3_ring];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}