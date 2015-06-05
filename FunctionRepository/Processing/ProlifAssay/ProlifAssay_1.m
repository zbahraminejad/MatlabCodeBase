cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\2013-04-18_DoseResponsePanel\';
timepoint = '24hr';

rows = [1:8];
cols = [1:12];
sites = [0:3];

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
cellcountmatrix=zeros(numrows,numcols,numsites);
timepath = [path,timepoint,'\'];

nucroot = '_ECFP_';
nucr=12;
minnucarea=pi*(nucr/4)^2;
frame=1;

for rowidx=1:numrows
    for colidx=1:numcols
        for siteidx=1:numsites
row=rows(rowidx); col=cols(colidx); site=sites(siteidx);
fprintf('row %0.0f, col %0.0f, site %0.0f\n',row,col,site);
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
DAs_or=single(imread([timepath,shot,nucroot,num2str(frame),'.tif']));

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
cellcountmatrix(rowidx,colidx,siteidx)=size(DAs_da,1);   %final readout
        end
    end
end
save([path,'cellcountmatrix_',timepoint,'_temp.mat'],'cellcountmatrix');
cd ..\Processing; %return to this directory