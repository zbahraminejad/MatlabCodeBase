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

SF=5;EF=50;
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
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    nuc_raw=log(single(imread([rawdir,nucname,num2str(frames(i)),'.tif'])));
    nuc_mask=blobdetector(nuc_raw,nucr);
    [~,numcells(i)]=bwlabel(nuc_mask);
    nuc_area(i)=median(cell2mat(struct2cell(regionprops(nuc_mask,'Area'))));
end
dims = size(nuc_mask);
height=dims(1); width=dims(2);
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
    nuc_raw=log(single(imread([rawdir,nucname,num2str(f),'.tif'])));
    YFP_raw=log(single(imread([rawdir,YFPname,num2str(f),'.tif'])));
    RFP_raw=log(single(imread([rawdir,RFPname,num2str(f),'.tif'])));
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector(nuc_raw,nucr);
    if i==firstgoodindex
        nuc_mask=segmentdeflections(nuc_mask,nucr,0.5); %0.5
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
    border=bwareaopen(nuc_mask,height*2+width*2);
    nuc_mask=logical(nuc_mask-border);
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,nuc_raw,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_bg=getbackground_H2B(nuc_mask,nuc_center,nuc_raw,nucr,0.25);
    nuc_density=nuc_density-nuc_bg;
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=nuc_density.*nuc_area;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==firstgoodindex
        daughters=ones(numcells,1)*NaN;
    else
        %%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mednuc=median(nuc_area);
        if numcells==0 || mednuc>blurthresh
            badframes(i)=1;
            imwrite(uint16(bwmorph(nuc_mask,'remove')),[maskdir,'\',nucedgename,num2str(f),'.tif']);
            continue;
        end
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        [jitx,jity]=getjitter_areamass(nuc_center,nuc_area,nuc_mass,tracedata{lastgoodframe}(:,1:4),dims);
        jitters(i,:)=jitters(lastgoodframe,:)+[jitx,jity];
        nuc_center(:,1)=nuc_center(:,1)+ones(numcells,1)*jitx;
        nuc_center(:,2)=nuc_center(:,2)+ones(numcells,1)*jity;
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        %[tracked,newdaughters,nuc_label]=adaptivetrack(tracedata{lastgoodframe}(:,1:3),initdata,nuc_raw,nuc_label,nucr);
        [tracked,newdaughters,nuc_label]=adaptivetrack_nextasterisk(tracedata{lastgoodframe}(:,1:4),initdata,nuc_raw,nuc_label,nucr,extractmask,jitx,jity);
        daughters=[daughters;newdaughters];
        numcells=size(tracked,1);
        nuc_center=tracked(:,[1 2]);
        nuc_area=tracked(:,3);
        nuc_mass=tracked(:,4);
        badframes(i)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    imwrite(uint16(extractmask),[maskdir,'\',nucedgename,num2str(f),'.tif']);
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(nuc_mass));
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    ring_label=getcytoring(nuc_label);
    ring_info=regionprops(ring_label,'PixelIdxList');
    nanvec=ones(numcells,1)*NaN; YFPnuc=nanvec; RFPnuc=nanvec; YFPring=nanvec; 
    for cc=cellid'
        YFPnuc(cc)=median(YFP_raw(nuc_info(cc).PixelIdxList));
        RFPnuc(cc)=median(RFP_raw(nuc_info(cc).PixelIdxList));
        ringall=YFP_raw(ring_info(cc).PixelIdxList);
         YFPring(cc)=mean(ringall(ringall>=median(ringall)));
    end
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [YFP_bg,RFP_bg]=getbackground_YR(nuc_mask,nuc_center,YFP_raw,RFP_raw,nucr,0.25);
    YFPnuc(cellid)=YFPnuc(cellid)-YFP_bg;
    RFPnuc(cellid)=RFPnuc(cellid)-RFP_bg;
    YFPring(cellid)=YFPring(cellid)-YFP_bg;
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracedata{i}=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,YFPnuc,RFPnuc,YFPring];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
%%% interpolate bad frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tracedata,jitters]=interpolateframes_adaptive(tracedata,jitters,badframes);
%%% compile tracestats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalcells=size(tracedata{end},1);
signal = zeros(totalcells,totalframes);
for i=1:totalframes
    %tempcell = find(~isnan(tracedata{i}(:,1)));
    %signal(tempcell,f) = bestsp{f}(tempcell,7);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}