function Immunostain_2(row,col,site)
row=5;col=1;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
%imagepath = 'E:\';
imagepath = 'H:\Images\';
shadingpath='H:\Images\ShadingImages\20140522 20xBin2\';

experimentpath='20131213 R-point CC\20140516 pRb Antibodies\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data\']);
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
name3='pRb';
name4='EdU';
nucedgename='nucedge';
moviebin=1;
if moviebin==1
    nucr=12; %MCF-10A:12 YT+:8
    debrisarea=100; %MCF-10A:200 YT+:300
    %nucr=20;
    %debrisarea=600;
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea; %MCF-10A:2000 YT+:2000
blobthreshold=-0.03; %MCF-10A:-0.03 YT+:-0.03
timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
% load([imagepath,experimentpath,'Raw\Serum\','DAPI_',num2str(site),'.mat']); sc1=shadingcorrection;
% load([imagepath,experimentpath,'Raw\Serum\','YFP_',num2str(site),'.mat']); sc2=shadingcorrection;
% load([imagepath,experimentpath,'Raw\Serum\','TxRed_',num2str(site),'.mat']); sc3=shadingcorrection;
% load([imagepath,experimentpath,'Raw\Serum\','Cy5_',num2str(site),'.mat']); sc4=shadingcorrection;
load([imagepath,experimentpath,'Raw\IBG\CFP_stain_',num2str(site),'.mat']); sc1=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\YFP_stain_',num2str(site),'.mat']); sc2=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\TexasRed_stain_',num2str(site),'.mat']); sc3=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\DAPI_stain_',num2str(site),'.mat']); sc4=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\Cy5_stain_',num2str(site),'.mat']); sc5=inferredbg;
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
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground=nuc_mask;
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %comment out for sparse YT analysis
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanmask=imdilate(foreground,strel('disk',nucr));
blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
real1=bgsubmasked_3(blur1,nanmask,1);
blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
real2=bgsubmasked_3(blur2,nanmask,1);
blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
real3=bgsubmasked_3(blur3,nanmask,1);
blur4=imfilter(raw4,fspecial('disk',3),'symmetric');
real4=bgsubmasked_3(blur4,nanmask,1);
%%% account for smearing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iqrmult=6;
% [real1,highfliers1]=bghighfliers(raw1,nanmask,nucr,iqrmult,10);
% [real2,highfliers2]=bghighfliers(raw2,nanmask,nucr,iqrmult,1);
% [real3,highfliers3]=bghighfliers(raw3,nanmask,nucr,iqrmult,10);
% nuc_label=bwlabel(nuc_mask);
% smearmask=highfliers1 | highfliers2 | highfliers3;
% smear_label=nuc_label.*smearmask; nuc_smeared=unique(smear_label(:)); nuc_smeared(1)=[];
% nuc_label(ismember(nuc_label,nuc_smeared))=0; nuc_mask=nuc_label>0;
%%% remove nucleoli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% feature_mask=smallblobdetector(log(raw2),12,-0.03,20,150);
% nuc_mask=logical(nuc_mask.*~feature_mask);
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=struct2cell(regionprops(nuc_mask,'Area','Centroid')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
if maskwrite
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label=bwlabel(nuc_mask);
numcells=numel(nuc_area);
nuc_info=regionprops(nuc_label,'PixelIdxList');
%ring_label=getcytoring(nuc_label);
ring_label=getcytoring_3(nuc_label,4,real2); %10xB1:4 20xB1:4
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig1=nanvec;
sig2=nanvec;
sig2ring_75th=nanvec; sig2ring_mode=nanvec;
sig3=nanvec;
sig4=nanvec;
for cc=1:numcells
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    sig4(cc)=mean(real4(nuc_info(cc).PixelIdxList));

    ringall=real2(ring_info(cc).PixelIdxList); 
    ringall(ringall>prctile(ringall,95))=[]; %hist(ringall,25);
    sig2ring_75th(cc)=prctile(ringall,75); %previously mean

    truering=ringall(ringall>25);

    if numel(truering)<50
         truering=ringall;
    end
    if numel(ringall)<2
         sig2ring_mode(cc)=NaN;
         continue;
    end             

    bmin=min(truering); bmax=max(truering); bstep=(bmax-bmin)/25; bins=bmin:bstep:bmax;
    [kval,xval]=ksdensity(truering,bins);
    maxidx=find(kval==max(kval),1); %first mode
    sig2ring_mode(cc)=xval(maxidx);
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
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig2ring_75th,sig2ring_mode,sig3,sig4];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig2ring_75th,sig2ring_50th,sig2ring_nobg,sig3,sig4];
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