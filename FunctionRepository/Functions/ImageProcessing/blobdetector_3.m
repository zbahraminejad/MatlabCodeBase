function nuc_mask = blobdetector_3(nuc_raw,nucr,threshold,debrisarea)

nuc_raw_sharp=imsharpenSean(nuc_raw,2,nucr*0.25); %added this on 2014/01/02
% threshold = -0.03;
% nuc_raw_sharp=nuc_raw;
sigma=0.5*nucr/sqrt(2); %default 0.75, bin1:0.75 bin2:0.5
h=sigma^2*fspecial('log',[round(nucr*2) round(nucr*2)],sigma); %laplacian of gaussian default window [nucr*2 nucr*2]
nuc_log=imfilter(nuc_raw_sharp,h,'symmetric');
nuc_mask=nuc_log<threshold; %higher picks up debris, lower misses parts of nuclei
nuc_mask=imfill(nuc_mask,'holes');
nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %bin1:4 bin2:2
nuc_mask=~bwmorph(~nuc_mask,'diag');
nuc_mask=~bwmorph(~nuc_mask,'bridge');
nuc_mask=bwareaopen(nuc_mask,debrisarea); 
% assessmask(nuc_raw,nuc_mask)

%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));  
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%