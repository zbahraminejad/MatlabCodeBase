function Immunostain_2(row,col,site)
%row=2;col=1;site=1;
%row=3;col=1;site=3;
row=3;col=1;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
%experimentpath='2013-11-21_IL2\Experiment_20131121\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131126\';
%experimentpath='2013-11-21_IL2\Experiment_20131204_pAKTpERK\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131209\';
%experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131213\';
experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131217\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data_backgroundtemplate_smear\']);
separatedirectories=0;
if separatedirectories
    rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir = [imagepath,experimentpath,'Mask\',shot];
else
    rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir = [imagepath,experimentpath,'Mask'];
end
bgdir=[imagepath,experimentpath,'Raw\'];
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='DAPI_'; %nuc
name2='FITC_';
name3='Cy3_';
name4='Cy5_';
nucedgename='nucedge_';
moviebin=1;
if moviebin==1
    nucr=12; %MCF-10A:12 YT+:8
    debrisarea=200; %MCF-10A:200 YT+:300
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
blobthreshold=-0.02; %YT+:-0.03
ignorenuclei=0;
timetotal=tic;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'stain.tif']));
raw2=single(imread([rawdir,name2,'stain.tif']));
raw3=single(imread([rawdir,name3,'stain.tif']));
raw4=single(imread([rawdir,name4,'stain.tif']));
%bgtemplate1=single(imread([bgdir,name1,'bgtemplate.tif']));
bgtemplate2=single(imread([bgdir,name2,'bgtemplate.tif']));
bgtemplate3=single(imread([bgdir,name3,'bgtemplate.tif']));
bgtemplate4=single(imread([bgdir,name4,'bgtemplate.tif']));
raw2=raw2-bgtemplate2; raw2(raw2<1)=1; raw2=log(imfilter(raw2,fspecial('disk',6),'symmetric'));
raw3=raw3-bgtemplate3; raw3(raw3<1)=1; raw3=log(imfilter(raw3,fspecial('disk',6),'symmetric'));
raw4=raw4-bgtemplate4; raw4(raw4<1)=1; raw4=log(imfilter(raw4,fspecial('disk',6),'symmetric'));
if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'stain.tif'];
end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
%boulderarea=2000; nuc_mask=blobdetector_excludelargeandwarped(raw1,nucr,blobthreshold,debrisarea,boulderarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask_withborder=nuc_mask;
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);

%%% detect smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smearrad=10;
smeararea=round(pi*6^2);
smearthreshold=-0.3;
cellmask=imdilate(nuc_mask_withborder,strel('disk',round(nucr*0.5)));
raw2withoutcells=raw2.*~cellmask; raw2withoutcells(cellmask)=median(raw2withoutcells(:));
raw3withoutcells=raw3.*~cellmask; raw3withoutcells(cellmask)=median(raw3withoutcells(:));
raw4withoutcells=raw4.*~cellmask; raw4withoutcells(cellmask)=median(raw4withoutcells(:));
raw2withoutcells=imresize(raw2withoutcells,0.1,'bicubic');
raw3withoutcells=imresize(raw3withoutcells,0.1,'bicubic');
raw4withoutcells=imresize(raw4withoutcells,0.1,'bicubic');
%smear2=smeardetector(raw2withoutcells,smearrad,smearthreshold,smeararea);
smear3=smeardetector(raw3withoutcells,smearrad,smearthreshold,smeararea);
%smear4=smeardetector(raw4withoutcells,smearrad,smearthreshold,smeararea);
%%% ignore nucleoli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ignorenuclei
    feature_mask=smallblobdetector(log(raw2),12,-0.03,20,150);
    nuc_mask_features=logical(nuc_mask.*~feature_mask);
else
    nuc_mask_features=nuc_mask;
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=struct2cell(regionprops(nuc_mask,raw1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_bg=centroidbackground(nuc_mask_withborder,nuc_center,raw1,nucr,3,0.25,50);
nuc_density=nuc_density-nuc_bg;
nuc_mass=nuc_density.*nuc_area;
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask_features,'remove');
if maskwrite
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=numel(nuc_mass);
nuc_info=regionprops(nuc_mask_features,'PixelIdxList');
%ring_label=getcytoring(nuc_label);
%ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig2=nanvec; sig2_mean=nanvec; sig2_max=nanvec;
sig3=nanvec; sig3_mean=nanvec; sig3_max=nanvec;
sig4=nanvec; sig4_mean=nanvec; sig4_max=nanvec;
%sig2ring=nanvec;
for cc=1:numcells
    sig2(cc)=median(raw2(nuc_info(cc).PixelIdxList));
    sig2_mean(cc)=mean(raw2(nuc_info(cc).PixelIdxList));
    sig2_max(cc)=max(raw2(nuc_info(cc).PixelIdxList));
    sig3(cc)=median(raw3(nuc_info(cc).PixelIdxList));
    sig3_mean(cc)=mean(raw3(nuc_info(cc).PixelIdxList));
    sig3_max(cc)=max(raw3(nuc_info(cc).PixelIdxList));
    sig4(cc)=median(raw4(nuc_info(cc).PixelIdxList));
    sig4_mean(cc)=mean(raw4(nuc_info(cc).PixelIdxList));
    sig4_max(cc)=max(raw4(nuc_info(cc).PixelIdxList));
    %ringall=YFP_raw(ring_info(cc).PixelIdxList);
     %sig2ring(cc)=mean(ringall(ringall>=median(ringall)));
end
%%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bg2=centroidbackground(nuc_mask_withborder,nuc_center,raw2,nucr,3,0.25,50); %AllRows:50
% bg3=centroidbackground(nuc_mask_withborder,nuc_center,raw3,nucr,3,0.25,50); %AllRows:50
% bg4=centroidbackground(nuc_mask_withborder,nuc_center,raw4,nucr,3,0.25,50); %AllRows:50
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
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig4,sig2_mean,sig3_mean,sig4_mean,sig2_max,sig3_max,sig4_max];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw2));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}