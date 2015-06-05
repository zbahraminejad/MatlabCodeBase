

% function VisualizeTracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%               VISUALIZE TRACKING              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Specify well/site position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row=3; col=6; site=1;
% shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%% Set Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projectpath='D:\Documents\Projects\';
% imagepath='D:\Images\';
% experimentpath='Michael\Roziglitazone\';
% separatedirectories=0;
% if separatedirectories==1
%     rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
%     %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
%     maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
% else
%     rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
%     maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
% end
% datadir=[projectpath,experimentpath,'Data\'];

%%%temp%%%%
% row=2;col=5 ;site=1;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projectpath='D:\Documents\Projects\';
% imagepath='G:\Projects\CellCycle_Differentiation\';
imagepath = 'G:\Zahra\';
%imagepath='E:\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
% experimentpath='11202014-Michael-CellCycle-48hour-glass\';
experimentpath='20150501-Zahra-Pulse\';%'12072014-Michael-48hr-CC\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
wellname = nameandsite(shot);
datadir=[imagepath,experimentpath,'Data\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    
    %rawdir=[imagepath,experi.mentpath,shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot];
    %maskdir=[imagepath,experimentpath,'Mask\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\'];
%     realdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracefile=[datadir,'tracedata_',shot,'_nolink','.mat'];
videofile=[datadir,'colortracking_nojitter_test_',shot,'.mat'];
load(tracefile,'tracedata','jitters','genealogy');
%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
showmask=1;
showcolor=1;
shownumber=0;
markerradius=4;
nucname='Cy5';
% nucname = [shot,'_CFP'];
% edgename='_nucedge';
edgename = [shot,'_nucedge'];
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
SF = 250;
EF = 300;
allframes = SF:EF;
numframes = length(allframes);
%%% Write to video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for f=1:numframes
for f=1:numframes
    fprintf('frame %0.0f\n',f);
    nucimage=imadjust(mat2gray(double(imread([rawdir,nucname,'_',num2str(allframes(f)),'.tif']))));
    if showmask
        edgeimage=logical(double(imread([maskdir,edgename,'_',num2str(allframes(f)),'.tif'])));
    else
        edgeimage=zeros(height,width);
    end
    centroidx=round(tracedata(:,allframes(f),1)-jitters(allframes(f),1));
    centroidy=round(tracedata(:,allframes(f),2)-jitters(allframes(f),2));
    centroidx(centroidx<1 | centroidx>width)=NaN;
    centroidy(centroidy<1 | centroidy>height)=NaN;
    centroidindex=sub2ind([height width],centroidy,centroidx);
    tracked=~isnan(centroidindex);
    %%% for cells being tracked, add a black marker at centroid
    if showcolor
        centroidmask=emptyframe;
        centroidmask(centroidindex(tracked))=1;
        centroidmask=imdilate(centroidmask,strel('disk',markerradius,0));
        nucimage(centroidmask==1)=0; %set to black
    end
    %%% detect mitosing cells
    trackbegins=firstframe==allframes(f);
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
    Rmother=emptyframe; Gmother=emptyframe; Bmother=emptyframe;
    if showcolor
        R(centroidindex(uniquenonmitosis))=cellcolor(uniquenonmitosis,1); R=imdilate(R,strel('disk',markerradius,0));
        G(centroidindex(uniquenonmitosis))=cellcolor(uniquenonmitosis,2); G=imdilate(G,strel('disk',markerradius,0));
        B(centroidindex(uniquenonmitosis))=cellcolor(uniquenonmitosis,3); B=imdilate(B,strel('disk',markerradius,0));
        %%% for mitosing cells in this frame, add mother's color
        Rmother(centroidindex(mitoses))=cellcolor(mitoses,1); Rmother=imdilate(Rmother,strel('disk',round(markerradius/2),0));
        Gmother(centroidindex(mitoses))=cellcolor(mitoses,2); Gmother=imdilate(Gmother,strel('disk',round(markerradius/2),0));
        Bmother(centroidindex(mitoses))=cellcolor(mitoses,3); Bmother=imdilate(Bmother,strel('disk',round(markerradius/2),0));
    end 
    %%% create RGB image
    imageframe=emptyRGB;
    imageframe(:,:,1)=nucimage+R+Rmother;
    imageframe(:,:,2)=nucimage+G+Gmother;
    imageframe(:,:,3)=nucimage+B+Bmother;
    imageframe(:,:,1)=imageframe(:,:,1).*~edgeimage;
    imageframe(:,:,2)=imageframe(:,:,2).*~edgeimage;
    imageframe(:,:,3)=imageframe(:,:,3).*~edgeimage;
    %%% add mask edge in green
    imageframe(:,:,2)=imageframe(:,:,2)+edgeimage;
    %%% overlay cell ID text
    if shownumber
        position=[centroidx centroidy];
        position=position(uniquecells,:);
        fontsize=12;
        textcolor='r';
        RGB=insertText(imageframe,position,find(uniquecells),'AnchorPoint','Center','FontSize',fontsize,'TextColor',textcolor,'BoxOpacity',0);
    else
        RGB=imageframe;
    end
    RGB(RGB>1)=1;
    RGB(RGB<0)=0;
    %%% calculate jitter calc
    frameminx=ceil(cropminx-jitters(f,1));
    framemaxx=floor(cropmaxx-jitters(f,1));
    frameminy=ceil(cropminy-jitters(f,2));
    framemaxy=floor(cropmaxy-jitters(f,2));
    frameheight=framemaxy-frameminy+1;
    framewidth=framemaxx-frameminx+1;
    %fprintf('cropheight=%0.0f\n',frameheight);
    %fprintf('cropwidth=%0.0f\n',framewidth);
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
    RGB=RGB(frameminy:framemaxy,frameminx:framemaxx,:);
    writeVideo(M,im2frame(RGB));
end
close(M);
% end