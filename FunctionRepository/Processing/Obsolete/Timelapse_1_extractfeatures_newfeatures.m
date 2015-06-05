%function Timelapse_1_extractfeatures(row,col,site)
%row='B';col='04';site='3';
row='D';col='05';site='4';
%row='C';col='05';site='1';

projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath = '2013-06-07_p21_cy2_deletions\Experiment_20130715\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
datadir = ([projectpath,experimentpath,'Data_Test\']);
rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
maskdir = [imagepath,experimentpath,'Mask_Test\',shot];
if ~exist(maskdir,'dir')
    mkdir(maskdir);
end

%SF=141;EF=160;
SF=1;EF=50;
nucr=12; %calculate avg nuclear radius with Measure_avgnucrad.m
nucname = 'CFP_'; %nuc
nucedgename = 'nucedge_';
YFPname = 'YFP_'; %DHB
RFPname = 'TexasRed_'; %p21

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=length(frames);
avgcellsize=pi*nucr^2;
debrissize=round(0.33*avgcellsize);
wellsss=cell(1,1,totalframes);  %change later to tracedata and cell(totalframes)
mednucsize=zeros(totalframes,1);
emptyframes=zeros(totalframes,1);
blurryframes=zeros(totalframes,1);

timetotal=tic;
for i=1:totalframes
    f=frames(i);
    fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_raw=log(single(imread([rawdir,nucname,num2str(f),'.tif'])));
    YFP_raw=log(single(imread([rawdir,YFPname,num2str(f),'.tif'])));
    RFP_raw=log(single(imread([rawdir,RFPname,num2str(f),'.tif'])));
    %%% extract image features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector(nuc_raw,nucr);
    %nuc_mask=segmentdeflections(nuc_mask,nucr,0.5); %third arg: min cell size ratio for calc
    nuc_mask=bwareaopen(nuc_mask,debrissize);
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=regionprops(nuc_mask,'Area','Centroid','PixelIdxList','Solidity','Eccentricity','MajorAxisLength');
    ring_label=getcytoring(nuc_label);
    ring_info=regionprops(ring_label,'PixelIdxList');
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_mask,'remove');
    imwrite(uint16(extractmask),[maskdir,'\',nucedgename,num2str(f),'.tif']);
    %%% check for empty frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numcells==0
        emptyframes(i)=1;
        continue;
    end
    %%% calculate background for each cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_bs,YFP_bs,RFP_bs]=getbackground(nuc_mask,nuc_info,nuc_raw,YFP_raw,RFP_raw,nucr,0.25);
    %%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tv=zeros(numcells,1);
    XX=tv; YY=tv; nuc_intensity=tv; nuc_area=tv; YFPnuc=tv; YFPring=tv; RFPnuc=tv;
    nuc_sol=tv; nuc_ecc=tv; nuc_mal=tv;
    for cc=1:numcells
        XX(cc)=nuc_info(cc).Centroid(1);  %x value of centroid
        YY(cc)=nuc_info(cc).Centroid(2);  %y value of centroid
        nuc_area(cc)=nuc_info(cc).Area;
        nuc_sol(cc)=nuc_info(cc).Solidity;
        nuc_ecc(cc)=nuc_info(cc).Eccentricity;
        nuc_mal(cc)=nuc_info(cc).MajorAxisLength;
        nuc_intensity(cc)=mean(nuc_raw(nuc_info(cc).PixelIdxList));
        RFPnuc(cc)=median(RFP_raw(nuc_info(cc).PixelIdxList));
        YFPnuc(cc)=median(YFP_raw(nuc_info(cc).PixelIdxList));
        allringpixels=YFP_raw(ring_info(cc).PixelIdxList);
            topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));
            YFPring(cc)=mean(topringpixels);
    end
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_intensity=nuc_intensity-nuc_bs;
    RFPnuc=RFPnuc-RFP_bs;
    YFPnuc=YFPnuc-YFP_bs;
    YFPring=YFPring-YFP_bs;
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wellsss{i}=[XX,YY,nuc_intensity,nuc_area,YFPnuc,RFPnuc,YFPring,nuc_sol,nuc_ecc,nuc_mal];
    mednucsize(i)=median(nuc_area);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
%%% detect bad frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blurthresh=median(mednucsize)*1.5;
blurryframes(mednucsize>blurthresh)=1;
badframes=emptyframes | blurryframes;
%%% calculate jitters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width] = size(nuc_mask);
[x,y]=calculatejitters_cellmatch(wellsss,badframes,height,width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'wellsss_',shot,'.mat'],'wellsss','badframes','x','y');

toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}