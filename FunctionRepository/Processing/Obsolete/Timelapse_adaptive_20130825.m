%function Timelapse_1_extractfeatures(row,col,site)
%row='B';col='04';site='3';
row='D';col='05';site='4';

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

SF=1;EF=50;
nucr=12; %calculate avg nuclear radius with Measure_avgnucrad.m
nucname = 'CFP_'; %nuc
nucedgename = 'nucedge_';
YFPname = 'YFP_'; %DHB
RFPname = 'TexasRed_'; %p21

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(totalframes,1)*NaN;
avgcellsize=pi*nucr^2;
debrissize=round(0.33*avgcellsize);
jitters=zeros(totalframes,2);
tracedata=cell(totalframes,1);  %change later to tracedata and cell(totalframes)
timetotal=tic;

%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(5,1); numcells=zeros(5,1);
for i=1:5
    CFP_raw=log(single(imread([rawdir,nucname,num2str(frames(i)),'.tif'])));
    nuc_mask=blobdetector(CFP_raw,nucr);
    [~,numcells(i)]=bwlabel(nuc_mask);
    nuc_area(i)=median(cell2mat(struct2cell(regionprops(nuc_mask,'Area'))));
end
dims = size(nuc_mask);
blurthresh=1.5*median(nuc_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=find(numcells>0 & nuc_area<blurthresh,1,'first');
if firstgoodindex>1
    for i=1:firstgoodindex-1
        badframes(i)=1;
        imwrite(uint16(zeros(dims)),[maskdir,'\',nucedgename,num2str(frames(i)),'.tif']);
    end
end
badframes(firstgoodindex)=0;
%%% derive cell trace signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CFP_raw=log(single(imread([rawdir,nucname,num2str(f),'.tif'])));
    YFP_raw=log(single(imread([rawdir,YFPname,num2str(f),'.tif'])));
    RFP_raw=log(single(imread([rawdir,RFPname,num2str(f),'.tif'])));
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector(CFP_raw,nucr);
    if i==firstgoodindex
        nuc_mask=segmentdeflections(nuc_mask,nucr,0.5);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_label,CFP_raw,'Area','Centroid','PixelIdxList')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    for cc=1:numcells
        CFPnuc=mean(CFP_raw(nuc_info(cc).PixelIdxList));
        YFPnuc=median(YFP_raw(nuc_info(cc).PixelIdxList));
        RFPnuc=median(RFP_raw(nuc_info(cc).PixelIdxList));
    end
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [CFP_bg,YFP_bg,RFP_bg]=getbackground(nuc_mask,nuc_center,CFP_raw,YFP_raw,RFP_raw,nucr,0.25);
    CFPnuc=CFPnuc-CFP_bg;
    YFPnuc=YFPnuc-YFP_bg;
    RFPnuc=RFPnuc-RFP_bg;
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=CFPnuc.*nuc_area;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>firstgoodindex
        %%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mednuc=median(nuc_area);
        if numcells==0 || mednuc>blurthresh
            badframes(i)=1;
            imwrite(uint16(bwmorph(nuc_mask,'remove')),[maskdir,'\',nucedgename,num2str(f),'.tif']);
            continue;
        end
        badframes(i)=0;
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jitxy=getjitter(nuc_center,nuc_mass,tracedata,dims);
        jitters(i,:)=jitters(find(badframes==0,1,'last'),:)+jitxy;
        nuc_center=nuc_center+ones(numcells,1)*jitxy;
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tracedata{i}=[nuc_center(:,1),nuc_center(:,2),nuc_mass,YFPnuc,RFPnuc,0];
        %%% track & correct merges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[nuc_label,tracedata]=adaptivetrack(nuc_label,tracedata);
        
        
        
        
        
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    imwrite(uint16(extractmask),[maskdir,'\',nucedgename,num2str(f),'.tif']);
    %%% add ring %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ring_label=getcytoring(nuc_label);
    ring_info=regionprops(ring_label,'PixelIdxList');
    YFPring=zeros(numcells,1); 
    for cc=1:numcells
        allringpixels=YFP_raw(ring_info(cc).PixelIdxList);
            topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));
            YFPring(cc)=mean(topringpixels);
    end
    YFPring=YFPring-YFP_bg;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
save([datadir,'tracedata_',shot,'.mat'],'tracedata','badframes','jitters');

toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}