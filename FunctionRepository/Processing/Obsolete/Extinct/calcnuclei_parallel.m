%% preparation
%rename files to be 001 002 003  [CHANGE]
close all;clear all; clear mex;

% if ~matlabpool('size')
%   matlabpool('open')
% end

time1=tic;
format short
ry1=0;ry2=0;rx1=0;rx2=0;  %to crop movie from start.  =0 if no cropping 
tempdir = 'h:\Documents\Timescape\20120807_12drugs';  %folder containing the movies folders [CHANGE]
rawdir = ([tempdir,'\Raw']);
%mkdir([tempdir,'\Data']);
SF=1;EF=208; %208
nucr=8; %use nucr=8 for MCF10A at 10x, nucr=9 for Hs68 and HeLa at 10x; nucr=13 for MCF10A at 20x

rows  = [3 4];
cols  = [1 2 3 4 5 6 7 8 9 10 11 12];
sites = [0 1];

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
shots=numrows*numcols*numsites;
parfor shot=1:shots
    siteidx=mod(shot,numsites);
    if siteidx==0
        siteidx=numsites;
    end
    site=sites(siteidx);
    colidx=mod(ceil(shot/numsites),numcols);
    if colidx==0
        colidx=numcols;
    end
    col=cols(colidx);
    rowidx=ceil(shot/(numcols*numsites));
    row=rows(rowidx);
    row,col,site
    main_body(row,col,site,SF,EF,tempdir,rawdir,nucr);
end
toc(time1)