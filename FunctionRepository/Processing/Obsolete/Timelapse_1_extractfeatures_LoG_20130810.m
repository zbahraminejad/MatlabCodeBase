%function Timelapse_1_extractfeatures_LoG(row,col,site)
%row='B';col='04';site='3';
row='D';col='05';site='4';
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
datadir = ([path,'Data\']);
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
rawdir = [path,'Raw\',shot,'\'];
cpdir = [path,'Processed\',shot];
if ~exist(cpdir,'dir')
    mkdir(cpdir);
end

SF=39;EF=230; %Sabrina20x:208 Steve20x:218 Steve10x:110 Steve&Sabrina:240
SF=229;
initF=SF;   %the first intended frame (only matters if I have to restart
%OldAxon: MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%IX-Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12;
nucname = 'CFP_'; %nuc
nucedgename = 'nucedge&ring_';
YFPname = 'YFP_'; %DHB
RFPname = 'TexasRed_'; %p21

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
continuation=0;
if continuation==1
    restartframe=76; %set to last frame successfully stored
    load([datadir,'wellsss_',shot,'_restart'],'wellsss');
    wellsssrestart=cell(1,1,EF-SF+1);
    wellsssrestart(1:restartframe-1)=wellsss(1:restartframe-1);
    wellsss=wellsssrestart; SF=restartframe;
else
    wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)  
end

timetotal=tic;
for f=SF:EF
    fprintf('frame %0.0f\n',f);
    timeframe=tic;

    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_raw=log(single(imread([rawdir,nucname,num2str(f),'.tif'])));
    YFP_raw=log(single(imread([rawdir,YFPname,num2str(f),'.tif'])));
    RFP_raw=log(single(imread([rawdir,RFPname,num2str(f),'.tif'])));
    %%% segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_pad=getnucmask_histanddeflection(nuc_bs,nucr,4,2); %previous segmentation
    nuc_mask=blobdetector_doublefilter(nuc_raw,nucr);
    %nuc_mask=detectdeflections(nuc_mask,nucr,1);
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=regionprops(nuc_mask,'Area','Centroid','PixelIdxList');
    ring_label=getcytoring(nuc_label);
    ring_info=regionprops(ring_label,'PixelIdxList');
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_mask,'remove');
    imwrite(uint16(extractmask),[cpdir,'\',nucedgename,num2str(f),'.tif']);
    %%% calculate background for each cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_bs,YFP_bs,RFP_bs]=getbackground(nuc_mask,nuc_info,nuc_raw,YFP_raw,RFP_raw,nucr);
    %%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tv=zeros(numcells,1);
    XX=tv; YY=tv; nuc_intensity=tv; nuc_area=tv; YFPnuc=tv; YFPring=tv; RFPnuc=tv;
    for cc=1:numcells
        XX(cc)=nuc_info(cc).Centroid(1);  %x value of centroid
        YY(cc)=nuc_info(cc).Centroid(2);  %y value of centroid
        nuc_area(cc)=nuc_info(cc).Area;
        nuc_intensity(cc)=median(nuc_raw(nuc_info(cc).PixelIdxList));
        RFPnuc(cc)=median(RFP_raw(nuc_info(cc).PixelIdxList));
        YFPnuc(cc)=median(YFP_raw(nuc_info(cc).PixelIdxList));
        allringpixels=YFP_raw(ring_info(cc).PixelIdxList);
            topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
            YFPring(cc)=mean(topringpixels);
    end
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_intensity=nuc_intensity-nuc_bs;
    RFPnuc=RFPnuc-RFP_bs;
    YFPnuc=YFPnuc-YFP_bs;
    YFPring=YFPring-YFP_bs;
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wellsss{f-initF+1}=[XX,YY,nuc_intensity,nuc_area,YFPnuc,RFPnuc,YFPring]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame
    save([datadir,'wellsss_', shot,'_test.mat'],'wellsss');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe)
end
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
cd([codepath,'Processing\']); %return to this directory