function nuc_mask = blobdetector(nuc_raw,nucr,threshold,debrisarea)
sigma=2*nucr/sqrt(2); %default 0.75
h=sigma^2*fspecial('log',[nucr*4 nucr*4],sigma); %laplacian of gaussian default window [nucr*2 nucr*2]
nuc_log=imfilter(nuc_raw,h,'symmetric');
nuc_mask=nuc_log<threshold; %higher picks up debris, lower misses parts of nuclei
nuc_mask=bwareaopen(nuc_mask,debrisarea); %bin1:0.75 bin2:0.5
nuc_mask=imdilate(nuc_mask,strel('disk',round(nucr)));


%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%