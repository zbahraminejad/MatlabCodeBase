function Timelapse_MSC_subBGfirst(row,col,site)
%row=6;col=7;site=1; %2 2 4
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
%shadingpath='H:\Images\ShadingImages\20140402_DAPI_CFP_YFP_TxRed_Cy5\';
%shadingpath='H:\Images\ShadingImages\20140410 DCYTC 20x\';
%shadingpath='H:\Images\ShadingImages\20140424 DCYTC 10x\';
shadingpath='H:\Images\ShadingImages\20140506 CYT 10x Long 200ul\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='20131213 R-point CC\20140407 2C CDK2 Hysteresis\';
%experimentpath='20131213 R-point CC\20140413 2C DHB Hysteresis\';
%experimentpath='Heewon\BreakBeforeMitosis\';
experimentpath='20131213 R-point CC\20140423 3C G1S 46i\';
%experimentpath='20131213 R-point CC\20140505 3C G1S TranxTrans\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=([projectpath,experimentpath,'Data_ksmode\']);
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
maskwrite=1;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
if separatedirectories==0
    maskdir=[maskdir,'\',shot,'_'];
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF=1;EF=103;
moviebin=1;
if moviebin==1
    nucr=12;
    debrisarea=100; %MCF-10A 10xBin1: H2B:100 NLS:200
    %nucr=20;
    %debrisarea=600;
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea; %MCF-10A: H2B:20 NLS:10
blobthreshold=-0.02;  %default 10xbin1=-0.02 20xbin2=-0.03;
%blobthreshold=-0.01; %NLS
jitterwindow=3*nucr; %default 10xbin1=48
ringcalc=1; %set to zero if ring is unneeded
name1='H2B_'; %nuc
name2='DHB_';
name3='Gem_';
namenucedge='nucedge_';
%adjustframe=0; %default
adjustframe=[41,42]; %G1S-46i
%adjustframe=[40,41,42,43]; %G1S-TranxTrans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(EF,1)*NaN;
if SF>1
    badframes(1:SF-1)=0;
end
jitters=zeros(EF,2);
blocksize=10000;
maxcellnum=blocksize;
tracedata=ones(maxcellnum,EF,12+ringcalc)*NaN; %default 7
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
timetotal=tic;
%[firstgoodindex,blurthreshhigh,numthresh,badframes,height,width]=timelapsesetup(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_3(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
%blurthreshlow=0;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
% load([shadingpath,'CFP','.mat'],'shadingcorrection'); sc1=shadingcorrection;
% load([shadingpath,'YFP','.mat'],'shadingcorrection'); sc2=shadingcorrection;
% load([shadingpath,'TxRed','.mat'],'shadingcorrection'); sc3=shadingcorrection;
load([imagepath,experimentpath,'Raw\','IBG_','H2B','.mat']); sc1=inferredbg;
load([imagepath,experimentpath,'Raw\','IBG_','DHB','.mat']); sc2=inferredbg;
load([imagepath,experimentpath,'Raw\','IBG_','Gem','.mat']); sc3=inferredbg;
% load([imagepath,experimentpath,'Raw\Serum\','DAPI_',num2str(site),'.mat']); sc1=shadingcorrection;
% load([imagepath,experimentpath,'Raw\Serum\','YFP_',num2str(site),'.mat']); sc2=shadingcorrection;
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw1=single(imread([rawdir,name1,num2str(f),'.tif'])); raw1=(raw1-bgcmos)./sc1;
    raw2=single(imread([rawdir,name2,num2str(f),'.tif'])); raw2=(raw2-bgcmos)./sc2;
    raw3=single(imread([rawdir,name3,num2str(f),'.tif'])); raw3=(raw3-bgcmos)./sc3;
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
    foreground=nuc_mask;
    if i==firstgoodindex
        nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %attempt on cells bigger than fraction (3rd arg)
        nuc_mask=excludelargeandwarped(nuc_mask,boulderarea); %default include
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
    border=bwareaopen(nuc_mask,height*2+width*2-4);
    nuc_mask=logical(nuc_mask-border);
    %%% background subtract: min %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    real1=bgsubmin(raw1,nucr);
    %real2=bgsubmin(raw2,nucr);
    %real3=bgsubmin(raw3,nucr);
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression=4;
    nanmask=imdilate(foreground,strel('disk',nucr/2));
    nanmaskcyto=imdilate(foreground,strel('disk',nucr));
    %real1=bgsubmasked_4(raw1,nanmask,1);
    blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
    real1=bgsubmasked_global(raw1,nanmask,1,compression);
    blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
    real1=bgsubmasked_global(raw1,nanmask,1,compression);
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if i>firstgoodindex && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)   
        fprintf('badframe: frame %0.0f\n',f);
        badframes(f)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
        end
        continue;
    end
    blurthreshhigh=1.1*mednuc;  %H2B:1.1 NLS:1.08
    blurthreshlow=0.8*mednuc;   %H2B:0.8 NLS:0.95
    numthresh=0.5*numcells;     %H2B:0.5 NLS:0.85
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
            %[reljitx,reljity]=getjitter_loners(nuc_center,relcenter,jitterwindow,[height width]);
            [reljitx,reljity,singlecount]=getjitter_loners_1(nuc_center,relcenter,jitterwindow,[height width]);
            if singlecount<30
                [reljitx,reljity]=detectGlobalShift(imfill(extractmask,'holes'),nuc_mask);
                reljitx=-reljitx;
                reljity=-reljity;
            end
        end
        jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        H2BvsNLS=1; %1:H2B 2:NLS
        debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_5(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),H2BvsNLS,debugpackage);
        badframes(f)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec; sig2median=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=mean(real1(nuc_info(cc).PixelIdxList)); %previously median
        sig2(n)=mean(real2(nuc_info(cc).PixelIdxList)); %previously mean
        sig2median(n)=median(real2(nuc_info(cc).PixelIdxList));
        %sig3(n)=mean(real3(nuc_info(cc).PixelIdxList)); %previously mean
        sig3(n)=median(real3(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        %ring_label=getcytoring(nuc_label);
        ring_label=getcytoring_3(nuc_label,4,real2); %10xB1:4 20xB1:4
        ring_info=regionprops(ring_label,'PixelIdxList');
        sig2ring=nanvec;
        sig2ring_75th=nanvec; sig2ring_50th=nanvec; sig2ring_nobg=nanvec; sig2ring_mode=nanvec;
        for n=1:numlivecells
            cc=cellid(n);
            ringall=real2(ring_info(cc).PixelIdxList);
            ringall(ringall>prctile(ringall,95))=[]; hist(ringall,25);
             sig2ring(n)=mean(ringall(ringall>=median(ringall))); %previously mean
             sig2ring_75th(n)=prctile(ringall,75); %previously mean
             sig2ring_50th(n)=median(ringall);
             
             truering=ringall(ringall>25);
             
             sig2ring_nobg(n)=median(truering); %20xBin1:75 20xBin2:50
                          
             if numel(truering)<50
                 %fprintf('truering low sample\n');
                 truering=ringall;
             end
             if numel(ringall)<2
                 sig2ring_mode(n)=NaN;
                 continue;
             end             
             
             bmin=min(truering); bmax=max(truering); bstep=(bmax-bmin)/25; bins=bmin:bstep:bmax;
             [kval,xval]=ksdensity(truering,bins);
             %[kval,xval]=hist(ringall,25);
             maxidx=find(kval==max(kval),1); %first mode
             sig2ring_mode(n)=xval(maxidx);
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc==1
        %tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1];
        %tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig2ring];
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2median,sig3,sig2ring,sig2ring_75th,sig2ring_50th,sig2ring_nobg,sig2ring_mode];
        %tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring];
        %tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,sig4];
        %tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig2ring];%Sabrina
    else
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
    end
    if maxcellnum-max(cellid)<blocksize
        tempdata=ones(blocksize,EF,12+ringcalc)*NaN; %default 7
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr,H2BvsNLS);
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
toc(timetotal);
clear all; clear mex;

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=imadjust(mat2gray(real2));
tempframe(:,:,3)=0;
tempframe(:,:,1)=extractmask;
figure,imshow(tempframe);
%}