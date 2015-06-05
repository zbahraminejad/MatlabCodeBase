%function moviemaker(row,col,site)
%%%%%%%% Make movie only for quick visual review %%%%%%%%%%%%%%%%%%%%%%
% Author: Mingyu Chung
% Last Revision: 10/19/2012 (Not verified)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all; clear mex;
row=5;col=12;site=1;
time1=tic;
format short
%path = 'h:\Documents\Timescape\20120807_12drugs\';  %folder containing the movies folders [CHANGE]
path = 'h:\Documents\MATLAB\2012-01-15_notreat\';
%datadir = ([path,'DataEval\']);
datadir = path;
rawdir = ([path,'picscropped\']);
SF=1;EF=120; %208
nucr=8; %use nucr=8 for MCF10A at 10x, nucr=9 for Hs68 and HeLa at 10x; nucr=13 for MCF10A at 20x


%% setup tempwell
movieName=[num2str(row),'_', num2str(col), '_', num2str(site)];  %[CHANGE]
disp(movieName);
M=avifile([datadir,movieName, 'PresentationMovie.avi'],'compression','none','fps',4);

for f=SF:EF
    time2=tic;

    %% reading images
    DAfile=[movieName,'_mRFP1_',num2str(f),'.tif'];
    YAfile=[movieName,'_EYFP_',num2str(f),'.tif'];
    DAs_or=single(imread([rawdir,DAfile]));
    YAs_or=single(imread([rawdir,YAfile]));  
    %% Make movie of cropped & processed images
    tempframe = imadjust(mat2gray(DAs_or));
    tempframe(:,:,2) = imadjust(mat2gray(YAs_or));
    tempframe(:,:,3) = zeros(size(YAs_or)); 
    M = addframe(M,im2frame(tempframe));

    toc(time2)
end
M=close(M);

size(DAs_or)  %display the size of the final cropped movie, to enter into Plot_Nuclei function
toc(time1)