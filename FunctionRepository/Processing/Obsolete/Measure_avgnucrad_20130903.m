function X_avgnucrad
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%path = 'h:\Documents\PCNA\20130516-H2B-PCNA-Wes\';
path = 'h:\Documents\CDK4\20130715\';
cpdir = [path,'CroppedProcessed\'];
nucroot = 'nuc_';
%IX_Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12; 
minnucarea=pi*(nucr/5)^2;
row='C';col='05';site='1';
frame=1;

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
%shot=[num2str(row),'_', num2str(col)];
shot=[row,'_',col,'_',site];
nuc_raw=single(imread([cpdir,shot,'\',nucroot,num2str(frame),'.tif']));
[height,width]=size(nuc_raw);

%%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_bs=bgsub(log(nuc_raw),10*nucr,0.05);

%%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=getnucmask_histsweep(nuc_bs,nucr);  %MC histogram sweep & concave detector
zeroborder=ones(height-2,width-2);
zeroborder=padarray(zeroborder,[1 1]);
nuc_mask=nuc_mask & zeroborder;   %necessary to imerode from edges
nuc_label=bwlabel(nuc_mask);
smoothmask=imopen(nuc_mask,strel('disk',round(nucr/5),0));
nuc_label=nuc_label.*smoothmask;  %maintains IDs, but smoothens & removes debris
nuc_info=regionprops(nuc_label,'Area','Centroid','PixelIdxList');

%%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(nuc_info,1);
nucexclude=zeros(numcells,1);
for k=1:numcells
    if nuc_info(k).Area < minnucarea
        nucexclude(k)=1;
    end
end
nucexclude=find(nucexclude);
nuc_info(nucexclude)=[];
numcells=size(nuc_info,1);

%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nuc_label(ismember(nuc_label,nucexclude))=0;
nucring=bwmorph(nuc_label,'remove');
tempfile=imadjust(mat2gray(nuc_bs));
tempfile(:,:,2)=nucring;
tempfile(:,:,3)=0;
imshow(tempfile);
%}
%%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AC=zeros(numcells,1);
for cc=1:numcells
    AC(cc)=nuc_info(cc).Area;
end

%%% calculate average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radii=sqrt(AC/pi);
figure;
hist(radii,100);
avgrad=round(mean(radii));
medrad=round(median(radii));
fprintf('average radius = %0.0f\n',avgrad);
fprintf('median radius = %0.0f\n',medrad);

cd ..\Processing; %return to this directory