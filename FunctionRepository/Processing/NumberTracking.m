 function ColorTracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%               VISUALIZE TRACKING              %%%%%%%%%%%%%%
%%%%%%%%%%%%%%                                               %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Specify well/site position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row=3; col=7; site=1;
%%% Set Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'E:\Michael\';
%imagepath='E:\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
% experimentpath='11202014-Michael-CellCycle-48hour-glass\';
experimentpath='20150126-gem-pparg-diff\';%'12072014-Michael-48hr-CC\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
wellname = nameandsite(shot);
datadir=[imagepath,experimentpath,'Data\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    
    %rawdir=[imagepath,experimentpath,shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot];
    %maskdir=[imagepath,experimentpath,'Mask\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\'];
%     realdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracefile=[datadir,'tracedata_',shot,'_withlinktest.mat'];
videofile=[datadir,'numbertracking_test_',shot];
load(tracefile,'tracedata','jitters','genealogy');
%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
markerradius=8;
nucname='mcherry';
dummyimage=double(imread([rawdir,nucname,'_1.tif']));
[height,width]=size(dummyimage);
imageframe=zeros(height,width,3);
emptyframe=zeros(height,width);
%%% Get first frame of each trace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numcells,numframes]=size(tracedata(:,:,1));
firstframe=ones(numcells,1)*NaN;
for c=1:numcells
    firstframe(c)=find(~isnan(tracedata(c,:,1)),1,'first');
end
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
    %%% note all cells being tracked
    tracked=~isnan(centroidindex);
    %%% ignore redundantly tracked cells
    sortedcentroidindex=sort(centroidindex);
    repeatedindex=find(diff(sortedcentroidindex)==0);
    redundantindex=sortedcentroidindex(repeatedindex);
    %centroidindex(ismember(centroidindex,redundantindex))=NaN;
    uniquecells=tracked & ~ismember(centroidindex,redundantindex);
    uniqueID=find(uniquecells);
    %centroidindex(mitoses)=NaN;
    %valid=~isnan(centroidindex);
    %%% write track ID on all unique cells
    imageframe(:,:,1)=nucimage;
    imageframe(:,:,2)=nucimage;
    imageframe(:,:,3)=nucimage;
    %posx=centroidx-markerradius; posy=centroidy-markerradius;
    %position=[posx posy];
    position=[centroidx centroidy];
    position=position(uniquecells,:);
    RGB=insertText(imageframe,position,uniqueID,'AnchorPoint','Center','FontSize',8,'TextColor','red','BoxOpacity',0);
    RGB(RGB>1)=1;
    writeVideo(M,im2frame(RGB));
end
close(M);
end