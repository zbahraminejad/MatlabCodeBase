function nuc_mask = blobdetector(nuc_raw,nucr,threshold)
sigma=0.5*nucr/sqrt(2); %default 0.75
h=sigma^2*fspecial('log',[nucr*1 nucr*1],sigma); %laplacian of gaussian
nuc_log=imfilter(nuc_raw,h,'symmetric');
nuc_mask=nuc_log>0.05; %higher picks up debris, lower misses parts of nuclei
%nuc_mask=imfill(nuc_mask,'holes');
%nuc_mask=imopen(nuc_mask,strel('disk',round(nucr/3),0));
avgcellsize=pi*nucr^2;
nuc_mask=~bwmorph(~nuc_mask,'diag');
nuc_mask=~bwmorph(~nuc_mask,'bridge');
nuc_mask=bwareaopen(nuc_mask,round(0.33*avgcellsize)); %default 0.33

%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%