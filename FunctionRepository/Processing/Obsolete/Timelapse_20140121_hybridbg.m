function Timelapse(row,col,site)
row='E';col='03';site='1';
%row=7;col=10;site=1;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
imagepath='H:\Images\';
%imagepath='E:\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130920\';
%experimentpath='2013-06-12_JY_BJ5_CycDFKBP_DHFRp21\';
%experimentpath='2012-12-11_Sabrina_MCF10A_DHFR-Chy-p21\';
%experimentpath='2013-02-22_Sabrina_MCF10Ap21null_DHFR-Chy-p21\';
%experimentpath='Kyuho\';
%experimentpath='HeeWon\';
%experimentpath='Sabrina\';
%experimentpath='Steve\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=([projectpath,experimentpath,'Data_20140121hybridbg\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',shot];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask'];
end
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
if separatedirectories==0
    maskdir=[maskdir,'\',shot,'_'];
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF=1;EF=230; %0715=230 0719=216 0905=220 0920=226 JY0612=208 SS1211=230 SS0222=149
moviebin=1;
if moviebin==1
    nucr=12;
    debrisarea=100; %MCF-10A: 200
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea;
blobthreshold=-0.02;  %default 10xbin1=-0.02
jitterwindow=3*nucr; %default 10xbin1=48
ringcalc=1; %set to zero if ring is unneeded
name1='CFP_'; %nuc
namenucedge='nucedge_';
name2='YFP_'; %DHB
name3='TexasRed_'; %p21
adjustframe=[37,211]; %20130715
%adjustframe=0; %20130719, SS1211, SS0222
%adjustframe=155; %20130831
%adjustframe=97; %20130905
%adjustframe=147; %20130920
%adjustframe=11; %Heewon 7_10_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(totalframes,1)*NaN;
jitters=zeros(totalframes,2);
blocksize=10000;
maxcellnum=blocksize;
tracedata=ones(maxcellnum,totalframes,6+ringcalc)*NaN;
%tracedata=ones(maxcellnum,totalframes,7+ringcalc)*NaN;
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
timetotal=tic;
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    raw1=log(single(imread([rawdir,name1,num2str(frames(i)),'.tif'])));
    nuc_mask=blobdetector(raw1,nucr,blobthreshold,debrisarea);
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
            imwrite(uint16(zeros(dims)),[maskdir,namenucedge,num2str(frames(i)),'.tif']);
        end
    end
end
badframes(firstgoodindex)=0;
blurthresh=1.08*nuc_area(firstgoodindex);
numthresh=0.5*numcells(firstgoodindex);
%%% derive cell trace signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw1=single(imread([rawdir,name1,num2str(f),'.tif']));
    raw2=single(imread([rawdir,name2,num2str(f),'.tif']));
    raw3=single(imread([rawdir,name3,num2str(f),'.tif']));
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
    if i==firstgoodindex
        nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %attempt on cells bigger than fraction (3rd arg)
        nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask_withborders=nuc_mask;
    nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
    border=bwareaopen(nuc_mask,height*2+width*2-4);
    nuc_mask=logical(nuc_mask-border);
    %%% background subtract %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nanmask=imdilate(nuc_mask_withborders,strel('disk',nucr));
    blocknum=10;
    raw2(nanmask)=NaN;
    raw3(nanmask)=NaN;
    bg2=blocksmooth(raw2,blocknum);
    bg3=blocksmooth(raw3,blocknum);
    real2=raw2-bg2;
    real3=raw3-bg3;
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,raw1,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if numcells<numthresh || (mednuc>blurthresh && i>firstgoodindex+1)
        fprintf('badframe: frame %0.0f\n',i);
        badframes(i)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[maskdir,namenucedge,num2str(f),'.tif']);
        end
        continue;
    end
    blurthresh=1.08*mednuc;  %cells can grow in size over time, so update each frame
    numthresh=0.5*numcells;
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
    %%% subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_bg=centroidbackground(nanmask,nuc_center,raw1,nucr,3,0.25,50);
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
        H2BvsNLS=1; %1:H2B 2:NLS
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_merge(lastgoodframe,f,tracedata,curdata,tracking,raw1,nuc_label,nanmask,nucr,jitters(i,:),H2BvsNLS);
        badframes(i)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; sig2=nanvec; sig3=nanvec; sig4=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
        sig3(n)=median(real3(nuc_info(cc).PixelIdxList));
        %sig4(n)=median(raw3(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        ring_label=getcytoring(nuc_label);
        ring_info=regionprops(ring_label,'PixelIdxList');
        sig2ring=nanvec; 
        for n=1:numlivecells
            cc=cellid(n);
            ringall=real2(ring_info(cc).PixelIdxList);
             sig2ring(n)=mean(ringall(ringall>=median(ringall)));
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc==1
        tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring];
        %tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,sig4];
        %tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig2ring];%Sabrina
    else
        tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3];
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
tracklengthcutoff=1; %default 5
excesstracks=(1:maxcellnum)'>maxcellid;
tracedata(excesstracks,:,:)=[];
tracking(excesstracks,:)=[];

%%% repair broken tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(tracedata,1);
emptytracks=zeros(numcells,1);
brokentracks=zeros(numcells,1);
tracestats=ones(numcells,2)*NaN;
for c=1:numcells
    if ~isempty(find(~isnan(tracedata(c,:,1)),1))
        tracestats(c,1)=find(~isnan(tracedata(c,:,1)),1,'first');
        tracestats(c,2)=find(~isnan(tracedata(c,:,1)),1,'last');        
    else
        emptytracks(c)=1;
    end
end
fin=find(tracestats(:,2)<totalframes-1);
nodaughters=zeros(numel(fin),1);
for i=1:numel(fin)
    nodaughters(i)=isempty(find(tracking(:,1)==fin(i),1));
end
fin=fin(logical(nodaughters));
finend=tracestats(fin,2);
[finend,tidx]=sort(finend,'descend');
fin=fin(tidx);
orphan=find(~isnan(tracestats(:,1)) & isnan(tracking(:,1)) & isnan(tracking(:,4)));
ostart=tracestats(orphan,1);
numorph=numel(orphan);
omass=ones(numorph,1)*NaN; ox=omass; oy=omass;
for i=1:numel(orphan)
    omass(i)=tracedata(orphan(i),ostart(i),4); ox(i)=tracedata(orphan(i),ostart(i),1); oy(i)=tracedata(orphan(i),ostart(i),2);
end
winrad=4*nucr;
for i=1:numel(fin)
    tmass=tracedata(fin(i),finend(i),4); tx=tracedata(fin(i),finend(i),1); ty=tracedata(fin(i),finend(i),2);
    c=find(ostart-finend(i)<4 & ostart-finend(i)>0 & abs(ox-tx)<winrad & abs(oy-ty)<winrad & abs((omass-tmass)/tmass)<0.3);
    if numel(c)>0
        c=c(find(ostart(c)==min(ostart(c)),1)); %multiple hits
        concatframes=ostart(c):totalframes;
        tracedata(fin(i),concatframes,:)=tracedata(orphan(c),concatframes,:);
        %tracedata(term(i),tend(i)+1,:)=(tracedata(term(i),tend(i),:)+tracedata(term(i),ostart(c),:))/2;
        framediff=ostart(c)-finend(i);
        if framediff>1
            datadiff=tracedata(fin(i),ostart(c),:) - tracedata(fin(i),finend(i),:);
            datastep=datadiff/framediff;
            for j=1:framediff-1
                tracedata(fin(i),finend(i)+j,:)=tracedata(fin(i),finend(i),:)+datastep*j;
            end
        end
        tracking((tracking(:,1)==orphan(c)),1)=fin(i); %reassign daughters
        ostart(c)=NaN; omass(c)=NaN; ox(c)=NaN; oy(c)=NaN;
        brokentracks(orphan(c))=1;
    end
end

%%% remove bad tracks (empty, broken, short, merged) %%%%%%%%%%%%%%%%%%%%%%
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
removetracks=find(emptytracks | brokentracks | shorttracks | mergetracks);
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
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=imadjust(mat2gray(raw1));
tempframe(:,:,3)=0;
tempframe(:,:,1)=extractmask;
figure,imshow(tempframe);

x6=tracedata(:,6,1);y6=tracedata(:,6,2);
x7=tracedata(:,7,1);y7=tracedata(:,7,2);

%}