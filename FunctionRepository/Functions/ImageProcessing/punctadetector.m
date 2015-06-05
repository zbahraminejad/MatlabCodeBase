function nuc_mask = punctadetector(nuc_raw,nucr,threshold,debrisarea,boulderarea)
%nuc_raw_sharp=imsharpenSean(nuc_raw,2,nucr*0.25); %added this on 2014/01/02
%nuc_raw_sharp=nuc_raw;
% sigma=0.75*nucr/sqrt(2); %default 0.75
% h=sigma^2*fspecial('log',[nucr*2 nucr*2],sigma); %laplacian of gaussian default window [nucr*2 nucr*2]
% nuc_log=imfilter(nuc_raw_sharp,h,'symmetric');
% nuc_mask=nuc_log<threshold; %higher picks up debris, lower misses parts of nuclei
% nuc_mask=imfill(nuc_mask,'holes');
% nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %bin1:4 bin2:2
% nuc_mask=~bwmorph(~nuc_mask,'diag');
% nuc_mask=~bwmorph(~nuc_mask,'bridge');

nuc_mask=nuc_raw>200;
%nuc_mask=bwareaopen(nuc_mask,debrisarea);
nuc_mask=imopen(nuc_mask,strel('disk',1,0));
antimask=imopen(nuc_mask,strel('disk',4,0));
nuc_mask=nuc_mask-antimask;

%antimask=bwareaopen(nuc_mask,boulderarea);
%nuc_mask=nuc_mask-antimask;

%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%dilated_mask=imdilate(nuc_mask,strel('disk',2,0));
%extractmask=bwmorph(dilated_mask,'remove');
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw_sharp));
tempframe(:,:,2)=extractmask;
%tempframe(:,:,2)=0;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%