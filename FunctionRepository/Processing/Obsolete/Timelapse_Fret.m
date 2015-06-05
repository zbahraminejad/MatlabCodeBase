function Timelapse_merge_Fret(row,col,site)
%row='A';col='02';site='1'; blur:12,50,63,69,79,96,126,127,132,144,192
%row='B';col='02';site='1';
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130920\';
%experimentpath='2013-06-12_JY_BJ5_CycDFKBP_DHFRp21\';
%experimentpath='2012-12-11_Sabrina_MCF10A_DHFR-Chy-p21\';
%experimentpath='2013-02-22_Sabrina_MCF10Ap21null_DHFR-Chy-p21\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=([projectpath,experimentpath,'Data_Merge\']);
rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
maskdir=[imagepath,experimentpath,'Mask\',shot];
maskwrite=1;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF=1;EF=220; %0715=230 0719=216 0905=220 0920=226 SS1211=230 SS0222=149
nucr=6; %BJ5-bin1=8
jitterwindow=6*nucr; %default 10xbin1=48
ringcalc=0; %set to zero if ring is unneeded
nucname='CFP_'; %nuc
nucedgename='nucedge_';
YFPname='YFP_'; %DHB
RFPname='TexasRed_'; %p21
blobthreshold=-0.02;  %default 10xbin1=-0.02
%adjustframe=[37,211]; %20130715
%adjustframe=0; %20130719
%adjustframe=155; %20130831
adjustframe=97; %20130905
%adjustframe=147; %20130920
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(totalframes,1)*NaN;
jitters=zeros(totalframes,2);
blocksize=10000;
maxcellnum=blocksize;
tracedata=ones(maxcellnum,totalframes,6+ringcalc)*NaN;
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
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
blurthresh=1.10*median(nuc_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=find(numcells>0 & nuc_area<blurthresh,1,'first');
if firstgoodindex>1
    for i=1:firstgoodindex-1
        badframes(i)=1;
        if maskwrite
            imwrite(uint16(zeros(dims)),[maskdir,'\',nucedgename,num2str(frames(i)),'.tif']);
        end
    end
end
badframes(firstgoodindex)=0;
blurthresh=1.10*nuc_area(firstgoodindex);
numthresh=0.1*numcells(firstgoodindex);
%%% derive cell trace signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_raw=single(imread([rawdir,nucname,num2str(f),'.tif']));
    YFP_raw=single(imread([rawdir,YFPname,num2str(f),'.tif']));
    RFP_raw=single(imread([rawdir,RFPname,num2str(f),'.tif']));
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector(log(nuc_raw),nucr,blobthreshold);
    if i==firstgoodindex
        nuc_mask=segmentdeflections(nuc_mask,nucr,0.5); %attempt on all non-debris
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
    if numcells<numthresh || mednuc>blurthresh
        badframes(i)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[maskdir,'\',nucedgename,num2str(f),'.tif']);
        end
        continue;
    end
    blurthresh=1.10*mednuc;  %cells can grow in size over time, so update each frame
    numthresh=0.1*numcells;
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_bg=getbackground_H2B(nuc_mask,nuc_center,nuc_raw,nucr,0.25);
    nuc_density=nuc_density-nuc_bg;
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=nuc_density.*nuc_area;
    curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>firstgoodindex
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        relprevx=tracedata(:,lastgoodframe,1)-jitters(lastgoodframe,1);
        relprevy=tracedata(:,lastgoodframe,2)-jitters(lastgoodframe,2);
        relcenter=[relprevx,relprevy];
        if ismember(f,adjustframe)
            [reljitx,reljity]=detectGlobalShift(imfill(extractmask,'holes'),nuc_mask);
            reljitx=-reljitx;
            reljity=-reljity;
        else
            [reljitx,reljity]=getjitter_loners(nuc_center,relcenter,jitterwindow,dims);
        end
        jitters(i,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        nuc_center(:,1)=nuc_center(:,1)+jitters(i,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(i,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_merge(lastgoodframe,f,tracedata,curdata,tracking,nuc_raw,nuc_label,nucr,jitters(i,:));
        badframes(i)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,'\',nucedgename,num2str(f),'.tif']);
    end
    %%% calculate FRET signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FRET_raw=YFP_raw;
    nuc_bg=bgsub(nuc_raw,10*nucr,0.05);
    FRET_bg=bgsub(FRET_raw,10*nucr,0.05);
    FRET_mask=blobdetector(log(FRET_raw),nucr,blobthreshold);
    [FRETjitx,FRETjity]=detectGlobalShift(imfill(extractmask,'holes'),imfill(FRET_mask,'holes'));
    FRET_bg=FRET_bg(1+FRETjitx:end-FRETjity,1+FRETjitx:end-FRETjitx);
    if FRETjitx>0
        padarray(FRET_bg,[0 FRETjitx],NaN,'post');
    else
        padarray(FRET_bg,[0 FRETjitx],NaN,'pre');
    end
    if FRETjity>0
        padarray(FRET_bg,[FRETjity 0],NaN,'post');
    else
        padarray(FRET_bg,[FRETjity 0],NaN,'pre');
    end
    FRET_ratio=FRET_bg./nuc_bg;
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; YFP_sig=nanvec; RFP_sig=nanvec; FRET_sig=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        YFP_sig(n)=median(YFP_raw(nuc_info(cc).PixelIdxList));
        RFP_sig(n)=median(RFP_raw(nuc_info(cc).PixelIdxList));
        FRET_sig(n)=median(FRET_sig(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        ring_label=getcytoring(nuc_label);
        ring_info=regionprops(ring_label,'PixelIdxList');
        YFP_sigring=nanvec; 
        for n=1:numlivecells
            cc=cellid(n);
            ringall=YFP_raw(ring_info(cc).PixelIdxList);
             YFP_sigring(n)=mean(ringall(ringall>=median(ringall)));
        end
    end
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc==1
        YFPprctile=10; RFPprctile=50;
    else
        YFPprctile=50; RFPprctile=50;
    end
    [YFP_bg,RFP_bg]=getbackground_12(nuc_mask,nuc_center,YFP_raw,RFP_raw,nucr,30,0.25,YFPprctile,RFPprctile);
    YFP_sig=YFP_sig-YFP_bg;
    RFP_sig=RFP_sig-RFP_bg;
    if ringcalc==1
        YFP_sigring=YFP_sigring-YFP_bg;
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc==1
        tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,YFP_sig,RFP_sig,YFP_sigring,FRET_sig];
    else
        tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,YFP_sig,RFP_sig,FRET_sig];
    end
    if maxcellnum-max(cellid)<3000
        tempdata=ones(blocksize,totalframes,6+ringcalc)*NaN;
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
%%% interpolate merged traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxcellid=max(cellid);
for i=1:maxcellid
    [validframes,numsets]=bwlabel(~isnan(tracedata(i,:,1)));
    if numsets>1
        for j=1:(numsets-1)
            mergepreframe=find(validframes==j,1,'last');
            mergepostframe=find(validframes==j+1,1,'first');
            mergeframes=(mergepreframe+1:mergepostframe-1)';
            datainit=squeeze(tracedata(i,mergepreframe,:));
            datadelta=squeeze(tracedata(i,mergepostframe,:)-tracedata(i,mergepreframe,:));
            totalstep=mergepostframe-mergepreframe;
            eachstep=mergeframes-mergepreframe;
            tracedata(i,mergeframes,:)=ones(numel(mergeframes),1)*datainit'+(eachstep*datadelta')/totalstep;
        end
    end
end
%%% interpolate bad frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tracedata,jitters]=interpolateframes_matrix(tracedata,jitters,badframes);
%%% remove excess tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracklengthcutoff=5;
excesstracks=(1:maxcellnum)'>maxcellid;
tracedata(excesstracks,:,:)=[];
tracking(excesstracks,:)=[];
%%% remove short tracks and merge tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracklength=ones(maxcellid,1)*NaN;
for c=1:maxcellid
    lastframe=find(~isnan(tracedata(c,:,1)),1,'last');
    firstframe=find(~isnan(tracedata(c,:,1)),1,'first');
    if isempty(lastframe)
        tracklength(c)=0;
    else
        tracklength(c)=lastframe-firstframe;
    end
end
mothers=unique(tracking(~isnan(tracking(:,1)),1));
shorttracks=isnan(tracking(:,1)) & ~ismember((1:maxcellid)',mothers) & tracklength<tracklengthcutoff;
mergetracks=~isnan(tracking(:,4));
removetracks=find(shorttracks | mergetracks);
tracedata(removetracks,:,:)=[];
tracking(removetracks,:)=[];
%%% update genealogy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genealogy=tracking(:,1);
for m=mothers'
    numremovedbefore=sum(removetracks<m);
    daughtersidx=genealogy==m;
    newm=m-numremovedbefore;
    genealogy(daughtersidx)=newm;
end

%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
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