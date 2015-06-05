function Timelapse_adaptive(row,col,site)
%row='G';col='08';site='2';

projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir = ([projectpath,experimentpath,'Data\']);
rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
maskdir = [imagepath,experimentpath,'Mask\',shot];
if ~exist(maskdir,'dir')
    mkdir(maskdir);
end

SF=1;EF=230;
nucr=7; %calculate avg nuclear radius with Measure_avgnucrad.m
nucname = 'CFP_'; %nuc
nucedgename = 'nucedge_';
YFPname = 'YFP_'; %DHB
RFPname = 'TexasRed_'; %p21

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(totalframes,1)*NaN;
jitters=zeros(totalframes,2);
tracedatacell=cell(totalframes,1);
blobthreshold=-0.02;  %default 10xbin1=-0.02. 10xbin2=-0.03. Try +/-0.01.
compress=1;     %retroactively convert bin1 to bin2
timetotal=tic;
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    nuc_raw=log(single(imread([rawdir,nucname,num2str(frames(i)),'.tif'])));
    nuc_mask=blobdetector(nuc_raw,nucr,blobthreshold);
    [~,numcells(i)]=bwlabel(nuc_mask);
    nuc_area(i)=median(cell2mat(struct2cell(regionprops(nuc_mask,'Area'))));
end
dims=size(nuc_mask);
height=dims(1); width=dims(2);
if compress
    height=round(height/2); width=round(width/2);
end
blurthresh=1.5*median(nuc_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=find(numcells>5 & nuc_area<blurthresh,1,'first');
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
    if compress
        nuc_raw=imresize(nuc_raw,0.5,'bilinear');
        YFP_raw=imresize(YFP_raw,0.5,'bilinear');
        RFP_raw=imresize(RFP_raw,0.5,'bilinear');
    end
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector(nuc_raw,nucr,blobthreshold);
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
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if numcells<=5 || mednuc>blurthresh
        badframes(i)=1;
        imwrite(uint16(bwmorph(nuc_mask,'remove')),[maskdir,'\',nucedgename,num2str(f),'.tif']);
        continue;
    end
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
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        relprevx=tracedatacell{lastgoodframe}(:,1)-jitters(lastgoodframe,1);
        relprevy=tracedatacell{lastgoodframe}(:,2)-jitters(lastgoodframe,2);
        relcenter=[relprevx,relprevy];
        [reljitx,reljity]=getjitter_loners(nuc_center,relcenter,nucr,dims);
        jitters(i,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        nuc_center(:,1)=nuc_center(:,1)+jitters(i,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(i,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        [tracked,newdaughters,nuc_label]=adaptivetrack(tracedatacell{lastgoodframe}(:,1:4),initdata,nuc_raw,nuc_label,nucr);
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
    tracedatacell{i}=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,YFPnuc,RFPnuc,YFPring];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
%%% interpolate bad frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badframes=find(badframes);
[tracedatacell,jitters]=interpolateframes_adaptive(tracedatacell,jitters,badframes);
%%% compile tracedata and tracestats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalcells=size(tracedatacell{end},1);
totalsignals=size(tracedatacell{end},2);
tracedata=ones(totalcells,totalframes,totalsignals)*NaN;
for f=1:totalframes
    presentcell = find(~isnan(tracedatacell{f}(:,1)));
    tracedata(presentcell,f,:) = tracedatacell{f}(presentcell,:);
end
tracestats=ones(totalcells,3)*NaN;
for c=1:totalcells
    goodframes=find(~isnan(tracedata(c,:,1)));
    if isempty(goodframes)  %un-merged cells
        continue;
    end
    firstframe=goodframes(1);
    lastframe=goodframes(end);
    tracestats(c,:)=[firstframe,daughters(c),lastframe];
end
%%% remove un-merged cells and short tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shortesttrack=5;
unmerged=find(isnan(tracestats(:,1)));
mothers=tracestats(:,2);
mothers=unique(mothers(~isnan(mothers)));
short=find(tracestats(:,3)-tracestats(:,1)<shortesttrack & isnan(tracestats(:,2)));
short=short(~ismember(short,mothers));
badtracks=[unmerged;short];
tracedata(badtracks,:,:)=[];
tracestats(badtracks,:)=[];
for m=mothers'
    numbadbefore=sum(badtracks<m);
    tempdaughters=find(tracestats(:,2)==m);
    newm=m-numbadbefore;
    tracestats(tempdaughters,2)=newm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
toc(timetotal);
clear all; clear mex;
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}