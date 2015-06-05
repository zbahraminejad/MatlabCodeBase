function Immunostain_2(row,col,site)
%row=1;col=5;site=3;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\Documents\Projects\';
%imagepath='D:\Images\';
imagepath='E:\';
%shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20131213 R-point CC\20140922 Fast46i p21cycA\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'Data20x\'];
biasdir=[imagepath,experimentpath,'Raw20x\Bias\'];
separatedirectories=1;
if separatedirectories
    rawdir = [imagepath,experimentpath,'Raw20x\',shot,'\'];
    maskdir = [imagepath,experimentpath,'Mask\',shot];
else
    rawdir = [imagepath,experimentpath,'Raw20x\',shot,'_'];
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
nucr=12;
debrisarea=200;
boulderarea=1200;
%blobthreshold=-0.03;
timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
[height,width]=size(bgcmos);
load([biasdir,name1,'_Stain_',num2str(site),'.mat']); bias1=bias;
load([biasdir,name2,'_Stain_',num2str(site),'.mat']); bias2=bias;
load([biasdir,name3,'_Stain_',num2str(site),'.mat']); bias3=bias;
load([biasdir,name4,'_Stain_',num2str(site),'.mat']); bias4=bias;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif'])); raw1=(raw1-bgcmos)./bias1;
raw2=single(imread([rawdir,name2,'_stain.tif'])); raw2=(raw2-bgcmos)./bias2;
raw3=single(imread([rawdir,name3,'_stain.tif'])); raw3=(raw3-bgcmos)./bias3;
raw4=single(imread([rawdir,name4,'_stain.tif'])); raw4=(raw4-bgcmos)./bias4;
if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%% standard segmentation
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground=nuc_mask;
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,0.5,debrisarea);
%}


%%% marker-based watershed: erosion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blurradius=3;
nuc_mask=threshmask(raw1,blurradius);
nuc_mask=markershed(nuc_mask,nucr*2/3);
foreground=nuc_mask;
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%anti_mask=bwareaopen(nuc_mask,boulderarea);
%nuc_mask=nuc_mask-anti_mask;

%{
%%% marker-based watershed: regionalmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=threshmask(raw1,3);
nuc_mask=markershed_regionalmax(nuc_mask,raw1,nucr);
foreground=nuc_mask;
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
%}

%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=imclearborder(nuc_mask);
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression=4;
nanmask=imdilate(foreground,strel('disk',nucr/2));
nanmaskcyto=imdilate(foreground,strel('disk',nucr));
blur1=imfilter(raw1,fspecial('gaussian',5),'symmetric');
real1=bgsubmasked_global(blur1,nanmask,1,compression);
blur2=imfilter(raw2,fspecial('gaussian',5),'symmetric');
real2=bgsubmasked_global(blur2,nanmaskcyto,1,compression);
blur3=imfilter(raw3,fspecial('gaussian',5),'symmetric');
real3=bgsubmasked_global(blur3,nanmask,1,compression);
blur4=imfilter(raw4,fspecial('gaussian',5),'symmetric');
real4=bgsubmasked_global(blur4,nanmaskcyto,1,compression);
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'bleedthroughrate_EdUtopRb.mat'],'bleedthroughrate');
real3bleedthrough=real4*bleedthroughrate(2)+bleedthroughrate(1);
real3=real3-real3bleedthrough;
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
%ring_label_org=getcytoring_3(nuc_label,4,real2); %10xB1:4 20xB1:4
%ring_info_org=regionprops(ring_label_org,'PixelIdxList');
innerrad=1; outerrad=5; %10xB1|20xB2: 1/5
ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig1=nanvec;
sig2=nanvec;
sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_fgmode=nanvec;
sig3=nanvec;
sig3ring_fgmedian=nanvec; sig3ring_fgmode=nanvec;
sig4=nanvec;
for cc=1:numcells
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    sig4(cc)=mean(real4(nuc_info(cc).PixelIdxList));

    ring2all=real2(ring_info(cc).PixelIdxList);
    ring3all=real3(ring_info(cc).PixelIdxList);
    %ringall(ringall>prctile(ringall,95))=[]; %hist(ringall,25);
    sig2ring_75th(cc)=prctile(ring2all,75);

    ring2foreground=ring2all(ring2all>0);
    ring3foreground=ring3all(ring2all>0);

    if numel(ring2foreground)<100
         ring2foreground=ring2all;
         ring3foreground=ring3all;
    end
    if numel(ring2all)>100
         sig2ring_fgmedian(cc)=nanmedian(ring2foreground);
         sig3ring_fgmedian(cc)=nanmedian(ring3foreground);
         numbins=25;
         sig2ring_fgmode(cc)=getmode(ring2foreground,numbins);
         sig3ring_fgmode(cc)=getmode(ring3foreground,numbins);
    end             
end
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode,sig3,sig3ring_fgmedian,sig3ring_fgmode,sig4];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}