function Immunostain_2(row,col,site)
row=2;col=8;site=6;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\Documents\Projects\';
imagepath='D:\Images\';
%imagepath='E:\';
%shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20140703 Cyclin Induction\20140923 CycEAind_EdU_pRb\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'Data\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
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
name3='cyclin';
name4='EdU';
secondname1='SecondHoechst';
secondname2='SecondpRb';

nucedgename='nucedge';
nucr=12;
debrisarea=200;
boulderarea=1200;
%nucr=20;
%debrisarea=700;

timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
[height,width]=size(bgcmos);
load([biasdir,name1,'_',num2str(site),'.mat']); bias1=bias;
load([biasdir,name2,'_',num2str(site),'.mat']); bias2=bias;
load([biasdir,name3,'_',num2str(site),'.mat']); bias3=bias;
load([biasdir,name4,'_',num2str(site),'.mat']); bias4=bias;
load([biasdir,secondname1,'_',num2str(site),'.mat']); biassec1=bias;
load([biasdir,secondname2,'_',num2str(site),'.mat']); biassec2=bias;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif'])); raw1=(raw1-bgcmos)./bias1;
raw2=single(imread([rawdir,name2,'_stain.tif'])); raw2=(raw2-bgcmos)./bias2;
raw3=single(imread([rawdir,name3,'_stain.tif'])); raw3=(raw3-bgcmos)./bias3;
raw4=single(imread([rawdir,name4,'_stain.tif'])); raw4=(raw4-bgcmos)./bias4;
rawsec1=single(imread([rawdir,secondname1,'_stain.tif'])); rawsec1=(rawsec1-bgcmos)./biassec1;
rawsec2=single(imread([rawdir,secondname2,'_stain.tif'])); rawsec2=(rawsec2-bgcmos)./biassec2;

if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% marker-based watershed: erosion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blurradius=3;
nuc_mask=threshmask(raw1,blurradius);
nuc_mask=markershed(nuc_mask,round(nucr*2/3));
foreground=nuc_mask;
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% segment second nuclear image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask2=threshmask(rawsec1,3);
nuc_mask2=markershed(nuc_mask2,nucr*2/3);
%%% align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx,reljity]=registerimages(nuc_mask,nuc_mask2);
reljitx=round(reljitx); reljity=round(reljity);
[nuc_mask,~]=cropboth(nuc_mask,nuc_mask2,reljitx,reljity);
[foreground,~]=cropboth(foreground,nuc_mask2,reljitx,reljity);
nuc_mask=logical(nuc_mask);
[raw1,~]=cropboth(raw1,rawsec1,reljitx,reljity);
[raw2,~]=cropboth(raw2,rawsec2,reljitx,reljity);
[raw3,~]=cropboth(raw3,rawsec2,reljitx,reljity);
[raw4,rawsec2]=cropboth(raw4,rawsec2,reljitx,reljity);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=imclearborder(nuc_mask);
%%% subtract EdU from pRb (both AF647) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blursec2=imfilter(rawsec2,fspecial('gaussian',5),'symmetric');
blur4=imfilter(raw4,fspecial('gaussian',5),'symmetric');
correctsec2=blursec2-blur4;
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression=4;
nanmask=imdilate(foreground,strel('disk',nucr/2));
nanmaskcyto=imdilate(foreground,strel('disk',nucr));
%blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
real1=bgsubmasked_global(raw1,nanmask,1,compression);
blur2=imfilter(raw2,fspecial('gaussian',5),'symmetric');
real2=bgsubmasked_global(blur2,nanmaskcyto,1,compression);
%blur3=imfilter(raw3,fspecial('gaussian',5),'symmetric');
real3=bgsubmasked_global(raw3,nanmaskcyto,1,compression);
real4=bgsubmasked_global(raw4,nanmask,1,compression);
realsec2=bgsubmasked_global(correctsec2,nanmask,1,compression);
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
innerrad=1; outerrad=5; %10xB1/20xB2: 1/5
ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig1=nanvec;
sig2=nanvec;
sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_fgmode=nanvec;
sig3=nanvec;
sig3ring_75th=nanvec; sig3ring_fgmedian=nanvec; sig3ring_fgmode=nanvec;
sig4=nanvec;
sigsec2=nanvec;
for cc=1:numcells
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=mean(real3(nuc_info(cc).PixelIdxList));
    sig4(cc)=mean(real4(nuc_info(cc).PixelIdxList));
    sigsec2(cc)=mean(realsec2(nuc_info(cc).PixelIdxList));

    ring2all=real2(ring_info(cc).PixelIdxList);
    ring3all=real3(ring_info(cc).PixelIdxList);
    sig2ring_75th(cc)=prctile(ring2all,75);
    sig3ring_75th(cc)=prctile(ring3all,75);

    ring2foreground=ring2all(ring2all>0);
    ring3foreground=ring3all(ring2all>0);

    if numel(ring2foreground)<50
         ring2foreground=ring2all;
         ring3foreground=ring3all;
    end
    if numel(ring2all)<50
         continue;
    end             
    sig2ring_fgmedian(cc)=nanmedian(ring2foreground);
    sig3ring_fgmedian(cc)=nanmedian(ring3foreground);
    numbins=25;
    sig2ring_fgmode(cc)=getmode(ring2foreground,numbins);
    sig3ring_fgmode(cc)=getmode(ring3foreground,numbins);
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_p21dCy1dKtop21.mat'],'bleedthroughrate');
% bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
% sig4=sig4-bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode,sig3,sig3ring_75th,sig3ring_fgmedian,sig3ring_fgmode,sig4,sigsec2];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
[height,width]=size(nuc_mask);
tempframe=zeros(height,width,3);
tempframe=imadjust(mat2gray(raw1));
%tempframe(:,:,1)=nuc_mask;
tempframe(:,:,2)=extractmask;
%tempframe(:,:,2)=nuc_mask2;
tempframe(:,:,3)=0;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}