function Immunostain_2(row,col,site)
%row=2;col=1;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
%experimentpath='2013-11-21_IL2\Experiment_20131121\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131126\';
%experimentpath='2013-11-21_IL2\Experiment_20131204_pAKTpERK\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131209\';
%experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131213\';
%experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131217\';
experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data_bgtemplateoffsetcalcblock\']);
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
name1='DAPI'; %nuc
name2='FITC';
name3='Cy3';
name4='Cy5';
nucedgename='nucedge';
moviebin=1;
if moviebin==1
    nucr=12; %MCF-10A:12 YT+:8
    debrisarea=200; %MCF-10A:200 YT+:300
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
blobthreshold=-0.03; %YT+:-0.03
removenuclei=0;
timetotal=tic;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif']));
raw2=single(imread([rawdir,name2,'_stain.tif']));
raw3=single(imread([rawdir,name3,'_stain.tif']));
raw4=single(imread([rawdir,name4,'_stain.tif']));

if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
%boulderarea=2000; nuc_mask=blobdetector_excludelargeandwarped(raw1,nucr,blobthreshold,debrisarea,boulderarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% calculate background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cytosol_mask=imdilate(nuc_mask,strel('disk',nucr*2));
blur1=imfilter(raw1,fspecial('disk',10),'symmetric');
blur2=imfilter(raw2,fspecial('disk',10),'symmetric');
blur3=imfilter(raw3,fspecial('disk',10),'symmetric');
blur4=imfilter(raw4,fspecial('disk',10),'symmetric');
iqrmult=3;
highfliers1=maskhighfliers(blur1,nuc_mask,iqrmult);
highfliers2=maskhighfliers(blur2,nuc_mask,iqrmult);
highfliers3=maskhighfliers(blur3,nuc_mask,iqrmult);
highfliers4=maskhighfliers(blur4,nuc_mask,iqrmult);
empty1=~cytosol_mask & ~highfliers1; empty1=empty1(:);
empty2=~cytosol_mask & ~highfliers2; empty2=empty2(:);
empty3=~cytosol_mask & ~highfliers3; empty3=empty3(:);
empty4=~cytosol_mask & ~highfliers4; empty4=empty4(:);
bgtemplate1=load([datadir,'Background_',name1,'.mat']);
bgtemplate2=load([datadir,'Background_',name2,'.mat']);
bgtemplate3=load([datadir,'Background_',name3,'.mat']);
bgtemplate4=load([datadir,'Background_',name4,'.mat']);
offset1=blur1(:)-bgtemplate1.meanrawrel(:,site);
offset2=blur2(:)-bgtemplate2.meanrawrel(:,site);
offset3=blur3(:)-bgtemplate3.meanrawrel(:,site);
offset4=blur4(:)-bgtemplate4.meanrawrel(:,site);
offset1=median(offset1(empty1));
offset2=median(offset2(empty2));
offset3=median(offset3(empty3));
offset4=median(offset4(empty4));
bg1=offset1+bgtemplate1.meanrawrel(:,site); bg1=vec2mat(bg1,width);
bg2=offset2+bgtemplate2.meanrawrel(:,site); bg2=vec2mat(bg2,width);
bg3=offset3+bgtemplate3.meanrawrel(:,site); bg3=vec2mat(bg3,width);
bg4=offset4+bgtemplate4.meanrawrel(:,site); bg4=vec2mat(bg4,width);
%%% remove nucleoli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% feature_mask=smallblobdetector(log(raw2),12,-0.03,20,150);
% nuc_mask=logical(nuc_mask.*~feature_mask);
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real1=raw1-bg1;
real2=raw2-bg2;
real3=raw3-bg3;
real4=raw4-bg4;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label=bwlabel(nuc_mask);
smearmask=highfliers1 | highfliers2 | highfliers3 | highfliers4;
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
%ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig1=nanvec; sig1_mean=nanvec; sig1_max=nanvec;
sig2=nanvec; sig2_mean=nanvec; sig2_max=nanvec;
sig3=nanvec; sig3_mean=nanvec; sig3_max=nanvec;
sig4=nanvec; sig4_mean=nanvec; sig4_max=nanvec;
%sig2ring=nanvec;
for cc=1:numcells
    sig1(cc)=median(real1(nuc_info(cc).PixelIdxList));
    sig1_mean(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig1_max(cc)=max(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig2_mean(cc)=mean(real2(nuc_info(cc).PixelIdxList));
    sig2_max(cc)=max(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    sig3_mean(cc)=mean(real3(nuc_info(cc).PixelIdxList));
    sig3_max(cc)=max(real3(nuc_info(cc).PixelIdxList));
    sig4(cc)=median(real4(nuc_info(cc).PixelIdxList));
    sig4_mean(cc)=mean(real4(nuc_info(cc).PixelIdxList));
    sig4_max(cc)=max(real4(nuc_info(cc).PixelIdxList));
    %ringall=YFP_raw(ring_info(cc).PixelIdxList);
     %sig2ring(cc)=mean(ringall(ringall>=median(ringall)));
end
%%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bg2=centroidbackground(nuc_mask_withborder,nuc_center,raw2,nucr,3,0.25,50); %AllRows:50
% bg3=centroidbackground(nuc_mask_withborder,nuc_center,raw3,nucr,3,0.25,50); %AllRows:50
%bg4=centroidbackground(nuc_mask_withborder,nuc_center,raw4,nucr,3,0.25,50); %AllRows:50
% sig2=sig2-bg2;
% sig2_mean=sig2_mean-bg2;
% sig2_max=sig2_max-bg2;
% sig3=sig3-bg3;
% sig3_mean=sig3_mean-bg3;
% sig3_max=sig3_max-bg3;
% sig4=sig4-bg4;
% sig4_mean=sig4_mean-bg4;
% sig4_max=sig4_max-bg4;
%sig2ring(cellid)=sig2ring(cellid)-bg2;
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
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2_mean,sig3_mean,sig2_max,sig3_max];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4,sig1_mean,sig2_mean,sig3_mean,sig4_mean,sig1_max,sig2_max,sig3_max,sig4_max];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(real2));
tempframe(:,:,2)=0;
tempframe(:,:,3)=0;
imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}