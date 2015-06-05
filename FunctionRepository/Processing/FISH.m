function FISH(row,col,site)
%row=7;col=6;site=2;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'D:\Documents\Projects\';
imagepath = 'D:\Images\';
shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='Heewon\20141007 cMyc FISH\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'Data\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=1;
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
name1='DAPI';
name2='mRNA';
name3='EdU';
nucedgename='nucedge';
nucr=12;
debrisarea=100;
secondmaskarea=200;
boulderarea=1200;
timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
[height,width]=size(bgcmos);
load([biasdir,name1,'.mat']); bias1=bias;
load([biasdir,name2,'.mat']); bias2=bias;
load([biasdir,name3,'.mat']); bias3=bias;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=double(imread([rawdir,name1,'_stain.tif'])); raw1=(raw1-bgcmos)./bias1; raw1(raw1<1)=1;
raw2=double(imread([rawdir,name2,'_stain.tif'])); raw2=(raw2-bgcmos)./bias2; raw2(raw2<1)=1;
raw3=double(imread([rawdir,name3,'_stain.tif'])); raw3=(raw3-bgcmos)./bias3; raw3(raw3<1)=1;
if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% marker-based watershed: erosion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blurradius=3;
nuc_mask=threshmask(raw1,blurradius);
nuc_mask=markershed(nuc_mask,nucr*2/3);
foreground=nuc_mask;
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=secondthresh(raw1,blurradius,nuc_mask,secondmaskarea);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
%nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% Note border/large/warped objects for later removal %%%%%%%%%%%%%%%%%%%%
antiborder_mask=imclearborder(nuc_mask);
border_mask=nuc_mask-antiborder_mask;
antilargewarped_mask=excludelargeandwarped(nuc_mask,boulderarea);
largewarped_mask=nuc_mask-antilargewarped_mask;
badobjects_mask=border_mask | largewarped_mask;
nuc_label=bwlabel(nuc_mask);
badIDs=unique(nuc_label.*badobjects_mask);
badIDs(1)=[];
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression=4;
nanmask=imdilate(foreground,strel('disk',nucr/2));
nanmaskcyto=imdilate(foreground,strel('disk',nucr));
real1=bgsubmasked_global(raw1,nanmask,1,compression);
real2global=bgsubmasked_global(raw2,nanmask,1,compression);
real3=bgsubmasked_global(raw3,nanmaskcyto,1,compression);
%%% segment puncta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real2=imtophat(raw2,strel('disk',2,0));
puncta_mask=real2>1000;
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
nuc_info=regionprops(nuc_label,'PixelIdxList');
thickenradius=2*nucr;
outer_label=labelthicken(nuc_label,thickenradius);
numcells=max(outer_label(:));
nuc_area=nuc_area(1:numcells);
nuc_center=nuc_center(1:numcells,:);
badIDs(badIDs>numcells)=[];
outer_info=regionprops(outer_label,'PixelIdxList');
puncta_label=outer_label.*puncta_mask;
puncta_info=regionprops(puncta_label,'PixelIdxList','Area');
punctacount=ones(numcells,1)*NaN;
for i=1:numcells
    tempmask=puncta_label==i;
    [~,tempcount]=bwlabel(tempmask,4);
    punctacount(i)=tempcount;
end
%%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanvec=ones(numcells,1)*NaN;
sig1=nanvec;
sig3=nanvec;
sig2totalglobal=nanvec;
sig2total=nanvec;
sig2puncint=nanvec;
sig2puncarea=nanvec;
sig2punccount=nanvec;
for cc=1:numcells
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig3(cc)=mean(real3(nuc_info(cc).PixelIdxList));
    sig2totalglobal(cc)=sum(real2global(outer_info(cc).PixelIdxList));
    sig2total(cc)=sum(real2(outer_info(cc).PixelIdxList));
    sig2puncint(cc)=sum(real2(puncta_info(cc).PixelIdxList));
    sig2puncarea(cc)=puncta_info(cc).Area;
    sig2punccount(cc)=punctacount(cc);
end
%%% remove bad object data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area(badIDs)=[];
nuc_center(badIDs,:)=[];
sig1(badIDs)=[];
sig3(badIDs)=[];
sig2totalglobal(badIDs)=[];
sig2total(badIDs)=[];
sig2puncint(badIDs)=[];
sig2puncarea(badIDs)=[];
sig2punccount(badIDs)=[];
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_p21dCy1dKtop21.mat'],'bleedthroughrate');
% bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
% sig4=sig4-bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig3,sig2totalglobal,sig2total,sig2puncint,sig2puncarea,sig2punccount];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
big_mask=imdilate(puncta_mask,strel('disk',1));
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}