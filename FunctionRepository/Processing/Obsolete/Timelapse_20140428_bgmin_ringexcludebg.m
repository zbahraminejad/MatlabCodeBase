function Timelapse_MSC_subBGfirst(row,col,site)
%row=3;col=7;site=1;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
%shadingpath='H:\Images\ShadingImages\20140402_DAPI_CFP_YFP_TxRed_Cy5\';
shadingpath='H:\Images\ShadingImages\20140410 DCYTC 20x\';
%shadingpath='H:\Images\ShadingImages\20140424 DCYTC 10x\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='20131213 R-point CC\20140407 2C CDK2 Hysteresis\';
experimentpath='20131213 R-point CC\20140413 2C DHB Hysteresis\';
%experimentpath='20131213 R-point CC\20140423 3C G1S 46i\';
%experimentpath='Heewon\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=([projectpath,experimentpath,'Data\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,shot,'\',shot,'_']; %Heewon
    maskdir=[imagepath,experimentpath,'Mask\',shot];
    %maskdir=[imagepath,experimentpath,'Mask\'];
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
SF=1;EF=51; %0715=230 0719=216 0905=220 0920=226 JY0612=208 SS1211=230 SS0222=149 20140208_p21dCy1dK=229
moviebin=1;
if moviebin==1
    %nucr=12;
    %debrisarea=100; %MCF-10A 10xBin1: H2B:100 NLS:200
    nucr=20;
    debrisarea=600;
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea; %MCF-10A: H2B:20 NLS:10
%blobthreshold=-0.02;  %default 10xbin1=-0.02
blobthreshold=-0.03; %20xbin1=-0.03
jitterwindow=3*nucr; %default 10xbin1=48
ringcalc=1; %set to zero if ring is unneeded
name1='H2B_'; %nuc
name2='DHB_'; %DHB
name3='Gem_'; %p21
namenucedge='nucedge_';
%adjustframe=0; %default
%adjustframe=[101,120,224,227]; %20140208_p21dCy1dK
%adjustframe=[37,211]; %20130715
%adjustframe=14; %20130719, SS1211, SS0222
%adjustframe=155; %20130831
%adjustframe=97; %20130905
%adjustframe=147; %20130920
%adjustframe=[11,12,41,42]; %Heewon 4_3_1
%adjustframe=[2,60,61,94]; %CyclingCDK4iWee1iCDC25i
%adjustframe=[25,26,30,31,43,44]; %20140407_2c_CDK2_Hysteresis
adjustframe=[2,24,25,26,34,35,36,38,39,40,47,48,49]; %DHB_Hysteresis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(totalframes,1)*NaN;
jitters=zeros(totalframes,2);
blocksize=10000;
maxcellnum=blocksize;
tracedata=ones(maxcellnum,totalframes,8+ringcalc)*NaN; %default 7
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
timetotal=tic;
[firstgoodindex,blurthreshhigh,numthresh,badframes,height,width]=timelapsesetup(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
%[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_1(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
blurthreshlow=0;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([shadingpath,'CFP','.mat'],'shadingcorrection'); sc1=shadingcorrection;
load([shadingpath,'YFP','.mat'],'shadingcorrection'); sc2=shadingcorrection;
% load([shadingpath,'BG_',num2str(site),'.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
% load([shadingpath,'CFP_',num2str(site),'.mat'],'shadingcorrection'); sc1=shadingcorrection;
% load([shadingpath,'YFP_',num2str(site),'.mat'],'shadingcorrection'); sc2=shadingcorrection;
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw1=single(imread([rawdir,name1,num2str(f),'.tif'])); raw1=(raw1-bgcmos)./sc1;
    raw2=single(imread([rawdir,name2,num2str(f),'.tif'])); raw2=(raw2-bgcmos)./sc2;
    %raw3=single(imread([rawdir,name3,num2str(f),'.tif']));
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
    if i==firstgoodindex
        nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %attempt on cells bigger than fraction (3rd arg)
        nuc_mask=excludelargeandwarped(nuc_mask,boulderarea); %default include
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
    border=bwareaopen(nuc_mask,height*2+width*2-4);
    nuc_mask=logical(nuc_mask-border);
    %%% background subtract %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blur1=imfilter(raw1,fspecial('disk',nucr),'symmetric');
    blur2=imfilter(raw2,fspecial('disk',nucr),'symmetric');
    c1=min(blur1(:));
    c2=min(blur2(:));
    bg1=ones(size(nuc_mask))*c1;
    bg2=ones(size(nuc_mask))*c2;
    real1=raw1-bg1;
    real2=raw2-bg2;
    %real3=raw3-bg3;
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    %if numcells<numthresh || (mednuc>blurthreshhigh && i>firstgoodindex+1)
    %if numcells<numthresh || ((mednuc>blurthreshhigh || mednuc<blurthreshlow) && i>firstgoodindex+1)
    if i>firstgoodindex && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)   
        fprintf('badframe: frame %0.0f\n',i);
        badframes(i)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[maskdir,namenucedge,num2str(f),'.tif']);
        end
        continue;
    end
    blurthreshhigh=1.1*mednuc;  %NLS:1.08
    blurthreshlow=0.8*mednuc;   %NLS:0.95
    numthresh=0.5*numcells;     %NLS:0.85
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
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
            [reljitx,reljity]=getjitter_loners(nuc_center,relcenter,jitterwindow,[height width]);
        end
        jitters(i,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        nuc_center(:,1)=nuc_center(:,1)+jitters(i,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(i,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        H2BvsNLS=1; %1:H2B 2:NLS
        debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_4(lastgoodframe,f,tracedata,curdata,tracking,raw1,nuc_label,nucr,jitters(i,:),H2BvsNLS,debugpackage);
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
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec; sig4=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=mean(real1(nuc_info(cc).PixelIdxList)); %previously median
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList)); %previously mean
        %sig3(n)=mean(real3(nuc_info(cc).PixelIdxList)); %previously mean
    end
    if ringcalc==1
        %ring_label=getcytoring(nuc_label);
        ring_label=getcytoring_2(nuc_label,4,raw2);
        ring_info=regionprops(ring_label,'PixelIdxList');
        sig2ring_75th=nanvec; sig2ring_50th=nanvec; sig2ring_nobg=nanvec;
        for n=1:numlivecells
            cc=cellid(n);
            ringall=real2(ring_info(cc).PixelIdxList); %hist(ringall,50);
             sig2ring_75th(n)=prctile(ringall,75); %previously mean
             sig2ring_50th(n)=median(ringall);
             sig2ring_nobg(n)=median(ringall(ringall>75));
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc==1
        %tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig2ring];
        tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_50th,sig2ring_nobg];
        %tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,sig4];
        %tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig2ring];%Sabrina
    else
        tracedata(cellid,i,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3];
    end
    if maxcellnum-max(cellid)<3000
        tempdata=ones(blocksize,totalframes,8+ringcalc)*NaN; %default 7
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
toc(timetotal);
clear all; clear mex;

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw2));
tempframe(:,:,2)=imadjust(mat2gray(raw1));
tempframe(:,:,3)=0;
tempframe(:,:,1)=extractmask;
figure,imshow(tempframe);
%}