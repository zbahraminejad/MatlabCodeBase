function nuc_mask = blobdetector(nuc_raw,nucr,threshold,smallthresh,bigthresh)
sigma=0.25*nucr/sqrt(2); %default 0.75
h=sigma^2*fspecial('log',[nucr*0.5 nucr*0.5],sigma); %laplacian of gaussian default window [nucr*2 nucr*2]
nuc_log=imfilter(nuc_raw,h,'symmetric');
nuc_mask=nuc_log<threshold; %higher picks up debris, lower misses parts of nuclei
nuc_mask=imfill(nuc_mask,'holes');
nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %bin1:4 bin2:2
nuc_mask=~bwmorph(~nuc_mask,'diag');
nuc_mask=~bwmorph(~nuc_mask,'bridge');
nuc_mask=bwareaopen(nuc_mask,smallthresh);
antimask=bwareaopen(nuc_mask,bigthresh);
nuc_mask=nuc_mask-antimask;


%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%