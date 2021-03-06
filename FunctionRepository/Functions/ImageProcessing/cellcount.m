function numcells = cellcount(row,col,site,timepoint)
%row=1;col=1;site=3;timepoint='24hr';
cd ..\Functions; %change directory for function calls
nucroot = '_ECFP_';
nucr=12;
minnucarea=pi*(nucr/4)^2;
frame=1;

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
DAs_or=single(imread([timepoint,'Raw\',shot,nucroot,num2str(frame),'.tif']));

%%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[height,width]=size(DAs_or);
%blocknum=3;
%blockheight=ceil(height/blocknum);
%blockwidth=ceil(width/blocknum);
%DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
DAs_bs=bgsub(log(DAs_or),10*nucr,0.05);

%%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector
zeroborder=ones(height-2,width-2);
zeroborder=padarray(zeroborder,[1 1]);
DAs_pad=DAs_pad & zeroborder;   %necessary to imerode from edges
realnuc_la=bwlabel(DAs_pad);
smoothmask=imopen(DAs_pad,strel('disk',round(nucr/6),0));
realnuc_la=realnuc_la.*smoothmask;  %maintains IDs, but smoothens & removes debris
DAs_da=regionprops(realnuc_la,'Area','Centroid','PixelIdxList');

%%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(DAs_da,1);
nucexclude=zeros(numcells,1);
for k=1:numcells
    if DAs_da(k).Area < minnucarea
        nucexclude(k)=1;
    end
end
nucexclude=find(nucexclude);
DAs_da(nucexclude)=[];
numcells=size(DAs_da,1);

%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
realnuc_la(ismember(realnuc_la,nucexclude))=0;
nucring=bwmorph(realnuc_la,'remove');
tempfile=imadjust(mat2gray(DAs_bs));
tempfile(:,:,2)=nucring;
tempfile(:,:,3)=0;
imshow(tempfile);
%}
cd ..\Processing; %return to this directory