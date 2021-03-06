% Image extraction
function fixedcell(row,col,site)

% row = 3;col = 2;site = 1;
tic
imagepath = 'D:\Michael\';
%imagepath='E:\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
experimentpath='12052014-CellCycle-24hr-48hr-bins\'; 
datadir = [imagepath,experimentpath,'Data\'];
rawdir = [imagepath, experimentpath,'Raw\'];
biasdir=[rawdir,'Bias\'];
%load the bias correction
% load([shadingpath,'\CameraNoise_2160_bin1.mat'],'BG'); bgcmos = BG;
load([shadingpath,'\cameranoise_4sites.mat'],'smoothnoise'); bgcmos = smoothnoise;
% define sites to visit

% channels = {'hoechst','bodipy','pparg','cebpb'};
name1 = 'hoechst';
name2 = 'pparg';
name3 = 'h2b';
% name4 = 'cebpb';

nucr = 12;
blobthreshold = -0.03;
debrisarea = 250;
boulderarea = 2000;

% load([biasdir,name1,'_Stain_',num2str(site),'.mat']);bias1=bias;
% load([biasdir,name2,'_Stain_',num2str(site),'.mat']);bias2=bias;
% load([biasdir,name3,'_Stain_',num2str(site),'.mat']);bias3=bias;
% load([biasdir,name4,'_Stain_',num2str(4),'.mat']);bias4=bias;
% day = 3;row = 5; col = 9; site = 1;

% rawdir = [imagepath,experimentpath,'Raw\','Day',num2str(day),'\'];

shot = [num2str(row),'_',num2str(col),'_',num2str(site)];

raw1 = double(imread([rawdir,shot,'_',name1,'.tif'])); 
raw2 = double(imread([rawdir,shot,'_',name2,'.tif'])); 
raw3 = double(imread([rawdir,shot,'_',name3,'.tif'])); 
% raw4 = double(imread([rawdir,shot,'_',name4,'.tif'])); 

tempmask = threshmask_fixed(raw1,3);
% tempmask2 = threshmask_fixed(raw2,3);
tempmask = imdilate(tempmask,strel('disk',nucr));
% tempmask3 = threshmask_fixed(raw3,3);
% tempmask4 = threshmask_fixed(raw4,3);
bias1 = illumbias(raw1,tempmask,bgcmos);
bias2 = illumbias(raw2,tempmask,bgcmos);
bias3 = illumbias(raw3,tempmask,bgcmos);
% bias4 = illumbias(raw4,tempmask4,bgcmos);

raw1=(raw1-bgcmos)./bias1; raw2=(raw2-bgcmos)./bias2; raw3=(raw3-bgcmos)./bias3; %raw4=(raw4-bgcmos)./bias4;

%high pass filter to to block the big things
%         tempraw = imfilter(raw1,fspecial('disk',3),'symmetric');
%         tempraw = imadjust(mat2gray(raw1));
%         tempthresh = graythresh(tempraw);
%         tempmask = im2bw(tempraw,tempthresh);
% %         tempmask = bwareaopen(tempmask,100);
%         antimask = bwareaopen(tempmask,15000);
% %         raw1 = imfilter(raw1,fspecial('disk',3));        
% %         raw1 = imtophat(raw1,strel('disk',nucr));
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground = nuc_mask;

%         
%         nucgray = imadjust(mat2gray(raw1));
%         thresh = graythresh(nucgray);
%         nucgray(antimask) = nan;
%         nuc_mask = im2bw(nucgray,thresh);
%         nuc_mask = threshmask(raw1,3);
%         foreground = nuc_mask;
%         nuc_mask = markershed(nuc_mask,6);
nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
% assessmask(raw1,nuc_mask)
nuc_mask = imclearborder(nuc_mask);
nuc_label = bwlabel(nuc_mask);

% dapicents = cell2mat(struct2cell(regionprops(nuc_label,'Centroid')));
% dapicents = reshape(dapicents,2,length(dapicents)/2);
% dapicents = dapicents';         


% blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
% bodipygray = imadjust(mat2gray(blur2));
% threshold = graythresh(bodipygray);
% bodipymask = im2bw(bodipygray,threshold);
% 
% bodipycents = cell2mat(struct2cell(regionprops(bodipymask,'Centroid')));
% bodipycents = reshape(bodipycents,2,length(bodipycents)/2);
% bodipycents = bodipycents';
% 
% numdapi = length(dapicents); 
% numbodipy = length(bodipycents);
% distmat = zeros(numdapi,numbodipy);
% for ii = 1:numdapi
% for jj = 1:numbodipy
% d = sqrt((bodipycents(jj,1)-dapicents(ii,1))^2+(bodipycents(jj,2)-dapicents(ii,2))^2);
% distmat(ii,jj) = d;
% end
% end
% [~,index] = sort(distmat);
% bodipylabel = index(1,:);
%establish nanmask with foreground and bodipy foreground
nanmask = imdilate(foreground,strel('disk',nucr*2));
% nanmaskbodipy = bodipymask;
compression = 4;
% background subtraction to get 'real' intensity values

blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
% blur4=imfilter(raw4,fspecial('disk',3),'symmetric');
real1=bgsubmasked_global_2(blur1,nanmask,1,compression,50);
real2=bgsubmasked_global_2(blur2,nanmask,1,compression,50); % nanmaskbodipy for the special case of bodipy, otherwise needs to be changed
real3=bgsubmasked_global_2(blur3,nanmask,1,compression,50);
% real4=bgsubmasked_global_2(blur4,nanmask,11,compression,50);

[nuc_label,numcells]=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
nuc_mass=nuc_density.*nuc_area;
nuc_info=regionprops(nuc_label,'PixelIdxList');
sig1 = zeros(numcells,1);
sig2 = zeros(numcells,1);
sig3 = zeros(numcells,1);
for cell=1:numcells
sig1(cell)=mean(real1(nuc_info(cell).PixelIdxList));
sig2(cell)=median(real2(nuc_info(cell).PixelIdxList));
sig3(cell)=median(real3(nuc_info(cell).PixelIdxList));
end



% get bodipy total intensity
% bodipyarea = cell2mat(struct2cell(regionprops(bodipymask,'Area')));
% bodipyint = cell2mat(struct2cell(regionprops(bodipymask,real2,'MeanIntensity')));
% bodipysum = bodipyarea.*bodipyint;
% [B,I] = sort(bodipylabel);
% bodipysum = bodipysum(I); % sorted based on labeled cell

%         finalbodipysum = zeros(length(dapicents),1);
%         finalbodipymedian = zeros(length(dapicents),1);

% sig2sum = zeros(length(dapicents),1);
% sig2median = zeros(length(dapicents),1);
% for ii = 1:max(B)
% sig2sum(ii) = sum(bodipysum(B==ii));
% sig2median(ii) = nanmedian(bodipysum(B==ii));
% end

fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2sum,sig2median,sig3,sig4];
% save([datadir,'fixdata_',shot,'_','Day',num2str(day),'.mat'],'fixdata');
% disp([shot,'_Day',num2str(day)])

save([datadir,'fixdata_Day4_',shot,'.mat'],'fixdata');
disp(shot)
toc
% end

