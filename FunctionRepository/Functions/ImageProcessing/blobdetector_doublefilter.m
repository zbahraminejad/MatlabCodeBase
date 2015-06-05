function logmaskfinal = blobdetector(nuc_raw,nucr)
sigma=0.5*nucr/sqrt(2);
h=sigma^2*fspecial('log',[nucr*5 nucr*5],sigma); %laplacian of gaussian
nuc_log=imfilter(nuc_raw,h,'symmetric');
logmask=nuc_log<0;
logmask=imopen(logmask,strel('disk',nucr/3,0));

sigma=0.1*nucr/sqrt(2);
h=sigma^2*fspecial('log',[nucr*5 nucr*5],sigma);
nuc_log=imfilter(nuc_raw,h,'symmetric');
logmaskfine=nuc_log<0;
logmaskfine=~bwmorph(~logmaskfine,'diag');      %break connections
logmaskfine=~bwmorph(~logmaskfine,'bridge');
logmaskfine=bwareaopen(logmaskfine,round(pi*(nucr/3)^2));

logmask_label=bwlabel(logmask);
overlap=logmask_label.*logmaskfine;
overlap=unique(overlap);
overlap(1)=[]; %remove 0
logmaskfinal=ismember(logmask_label,overlap);
logmaskfinal=imfill(logmaskfinal,'holes');
logmaskfinal=bwareaopen(logmaskfinal,round(pi*(nucr*4/7)^2));

%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(logmaskfinal,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%