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
videofile=[datadir,'colortracking_nojitter_',shot,'.mat'];
load(tracefile,'tracedata','jitters','genealogy');
%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markerradius=8;
nucname='CFP';
dummyimage=double(imread([rawdir,nucname,'_1.tif']));
[height,width]=size(dummyimage);
emptyframe=zeros(height,width);
emptyRGB=zeros(height,width,3);
%%% Calculate max cropping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jitminx=min(jitters(:,1));
jitmaxx=max(jitters(:,1));
jitminy=min(jitters(:,2));
jitmaxy=max(jitters(:,2));
cropminx=1; cropmaxx=width; cropminy=1; cropmaxy=height;
if jitmaxx>0
    cropminx=1+jitmaxx;
end
if jitminx<0
    cropmaxx=width-abs(jitminx);
end
cropwidth=floor(cropmaxx-cropminx)+1;
if jitmaxy>0
    cropminy=1+jitmaxy;
end
if jitminy<0
    cropmaxy=height-abs(jitminy);
end
cropheight=floor(cropmaxy-cropminy)+1;
%%% Get first frame of each trace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numcells,numframes]=size(tracedata(:,:,1));
firstframe=ones(numcells,1)*NaN;
for c=1:numcells
    firstframe(c)=find(~isnan(tracedata(c,:,1)),1,'first');
end
%%% Assign colors to cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    centroidx(centroidx<1 | centroidx>width)=NaN;
    centroidy(centroidy<1 | centroidy>height)=NaN;
    centroidindex=sub2ind([height width],centroidy,centroidx);
    %%% for cells being tracked, add a black marker at centroid
    tracked=~isnan(centroidindex);
    centroidmask=emptyframe;
    centroidmask(centroidindex(tracked))=1;
    centroidmask=imdilate(centroidmask,strel('disk',markerradius,0));
    nucimage(centroidmask==1)=0; %set to black
    %%% detect mitosing cells
    trackbegins=firstframe==f;
    mitoses=tracked & trackbegins & ~isnan(genealogy);
    mothers=genealogy(mitoses);
    %%% for cells that are redundant or mitosing, leave black (don't color)
    sortedcentroidindex=sort(centroidindex);
    repeatedindex=find(diff(sortedcentroidindex)==0);
    redundantindex=sortedcentroidindex(repeatedindex);
    %centroidindex(ismember(centroidindex,redundantindex))=NaN;
    uniquecells=tracked & ~ismember(centroidindex,redundantindex);
    uniquenonmitosis=uniquecells & ~mitoses;
    %centroidindex(mitoses)=NaN;
    %valid=~isnan(centroidindex);
    %%% color all unique, non-mitosing cells by track ID color
    R=emptyframe; G=emptyframe; B=emptyframe;
    R(centroidindex(uniquenonmitosis))=cellcolor(uniquenonmitosis,1); R=imdilate(R,strel('disk',markerradius,0));
    G(centroidindex(uniquenonmitosis))=cellcolor(uniquenonmitosis,2); G=imdilate(G,strel('disk',markerradius,0));
    B(centroidindex(uniquenonmitosis))=cellcolor(uniquenonmitosis,3); B=imdilate(B,strel('disk',markerradius,0));
    %%% for mitosing cells in this frame, add mother's color
    Rmother=emptyframe; Gmother=emptyframe; Bmother=emptyframe;
    Rmother(centroidindex(mitoses))=cellcolor(mothers,1); Rmother=imdilate(Rmother,strel('disk',round(markerradius/2),0));
    Gmother(centroidindex(mitoses))=cellcolor(mothers,2); Gmother=imdilate(Gmother,strel('disk',round(markerradius/2),0));
    Bmother(centroidindex(mitoses))=cellcolor(mothers,3); Bmother=imdilate(Bmother,strel('disk',round(markerradius/2),0));
    %%% create RGB image
    imageframe=emptyRGB;
    imageframe(:,:,1)=nucimage+R+Rmother;
    imageframe(:,:,2)=nucimage+G+Gmother;
    imageframe(:,:,3)=nucimage+B+Bmother;
    imageframe(imageframe>1)=1;
    imageframe(imageframe<0)=0;
    %%% calculate jitter calc
    frameminx=ceil(cropminx-jitters(f,1));
    framemaxx=floor(cropmaxx-jitters(f,1));
    frameminy=ceil(cropminy-jitters(f,2));
    framemaxy=floor(cropmaxy-jitters(f,2));
    frameheight=framemaxy-frameminy+1;
    framewidth=framemaxx-frameminx+1;
    fprintf('cropheight=%0.0f\n',frameheight);
    fprintf('cropwidth=%0.0f\n',framewidth);
    if framewidth<cropwidth
        if framemaxx<width
            framemaxx=framemaxx+1;
        else
            frameminx=frameminx-1;
        end
    end
    if frameheight<cropheight
        if framemaxy<height
            framemaxy=framemaxy+1;
        else
            frameminy=frameminy-1;
        end
    end
    %%% crop frames to remove jitters
    imageframe=imageframe(frameminy:framemaxy,frameminx:framemaxx,:);
    writeVideo(M,im2frame(imageframe));
end
close(M);
end