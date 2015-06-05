function moviemaker(row,col,site)
%%%%%%%% Make movie only for quick visual review %%%%%%%%%%%%%%%%%%%%%%
% Author: Mingyu Chung
% Last Revision: 10/19/2012 (Not verified)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all; clear mex;

time1=tic;
format short
ry1=0;ry2=0;rx1=0;rx2=0;  %to crop movie from start.  =0 if no cropping 
path = 'h:\Documents\Timescape\20120807_12drugs\';  %folder containing the movies folders [CHANGE]
datadir = ([path,'DataEval\']);
rawdir = ([path,'Raw\']);
SF=1;EF=30; %208
nucr=8; %use nucr=8 for MCF10A at 10x, nucr=9 for Hs68 and HeLa at 10x; nucr=13 for MCF10A at 20x


%% setup tempwell
movieName=[num2str(row),'_', num2str(col), '_', num2str(site)];  %[CHANGE]
disp(movieName);
M=avifile([datadir,movieName, 'H2B.avi'],'compression','none','fps',4);
%% correcting jitters
x=zeros(1,EF); y=zeros(1,EF);

for f=SF:EF
    disp(f) %print to the screen the frame you are on
    time2=tic;

    %% reading images
    DAs_or=single(imread([rawdir,movieName, '_mRFP1_', num2str(f), '.tif']));%DA = dapi, reads image. 2D matrix of intensities for that frame

    %% correct multiple jitters
    %DI1=imadjust(mat2gray(log(DAs_or(1:925, :))));  %only use the first 925 rows of pixels bc sometimes images have black shadows/lines at the bottom
    DI1=imadjust(mat2gray(log(DAs_or)));
    if f>SF                   
        [xx,yy]=CalcJitter(DI0,DI1);
        x(f)=x(f-1)+xx;  
        y(f)=y(f-1)+yy;
        disp([x(f), y(f)]);
    end
    DI0=DI1;
    toc(time2);
end

cropamountx=ceil(max(abs(x))+1);  %use maximum jitter to find amount to crop  
cropamounty=ceil(max(abs(y))+1);  %use maximum jitter to find amount to crop         

for f=SF:EF
    time2=tic;

    %% reading images
    DAfile=[movieName,'_mRFP1_',num2str(f),'.tif'];
    DAs_or=single(imread([rawdir,DAfile]));

    %% cropping.  
    DAs_or=CropJitter(DAs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));  %75 is my guess at the biggest jitter in pixels in whole movie

    %% Make movie of cropped & processed images
    viewframe = imadjust(mat2gray(DAs_or));
    tempframe = viewframe;  %red; to make a movie incluing nuclear marker
    tempframe(:,:,2) = viewframe;   %green
    tempframe(:,:,3) = viewframe; %blue.  
    M = addframe(M,im2frame(tempframe));

    toc(time2)
end
M=close(M);

size(DAs_or)  %display the size of the final cropped movie, to enter into Plot_Nuclei function
toc(time1)
end