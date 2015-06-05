function ColorTracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%               VISUALIZE TRACKING              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Specify well/site position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row=2; col=3; site=1;
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%% Set Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\Documents\Projects\';
imagepath='D:\Images\';
experimentpath='Michael\OP9\';
rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
datadir=[projectpath,experimentpath,'Data\'];
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracefile=[datadir,'tracedata_',shot,'.mat'];
videofile=[datadir,'colortracking_',shot,'.mat'];
load(tracefile,'tracedata','jitters');
%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markerradius=6;
nucname='CFP';
dummyimage=double(imread([rawdir,nucname,'_1.tif']));
[height,width]=size(dummyimage);
imageframe=zeros(height,width,3);
emptyframe=zeros(height,width);
%%% Assign colors to cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numcells,numframes]=size(tracedata(:,:,1));
colormapsize=64;
colormap=jet(colormapsize); close(gcf);
colormapidx=ceil(rand(numcells,1)*colormapsize);
cellcolor=colormap(colormapidx,:);
%%% Initiate video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=VideoWriter(videofile,'Uncompressed AVI');
M.FrameRate=4;
open(M);
%%% Write to video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f=1:numframes
    fprintf('frame %0.0f\n',f);
    nucimage=imadjust(mat2gray(double(imread([rawdir,nucname,'_',num2str(f),'.tif']))));
    centroidx=round(tracedata(:,f,1)-jitters(f,1));
    centroidy=round(tracedata(:,f,2)-jitters(f,2));
    centroidindex=sub2ind([height width],centroidy,centroidx);
    %%% for cells being tracked, add a black marker at centroid
    valid=~isnan(centroidindex);
    centroidmask=emptyframe;
    centroidmask(centroidindex(valid))=1;
    centroidmask=imdilate(centroidmask,strel('disk',markerradius,0));
    nucimage(centroidmask==1)=0; %set to black
    %%% for cells that are redundantly tracked, leave black (don't color)
    sortedcentroidindex=sort(centroidindex);
    repeatedindex=find(diff(sortedcentroidindex)==0);
    redundantindex=sortedcentroidindex(repeatedindex);
    centroidindex(ismember(centroidindex,redundantindex))=NaN;
    valid=~isnan(centroidindex);
    %%% color all non-redundantly tracked cells by the 'cellcolor' legend
    R=emptyframe; G=emptyframe; B=emptyframe;
    R(centroidindex(valid))=cellcolor(valid,1); R=imdilate(R,strel('disk',markerradius,0));
    G(centroidindex(valid))=cellcolor(valid,2); G=imdilate(G,strel('disk',markerradius,0));
    B(centroidindex(valid))=cellcolor(valid,3); B=imdilate(B,strel('disk',markerradius,0));
    imageframe(:,:,1)=nucimage+R;
    imageframe(:,:,2)=nucimage+G;
    imageframe(:,:,3)=nucimage+B;
    imageframe(imageframe>1)=1;
    writeVideo(M,im2frame(imageframe));
end
close(M);
end