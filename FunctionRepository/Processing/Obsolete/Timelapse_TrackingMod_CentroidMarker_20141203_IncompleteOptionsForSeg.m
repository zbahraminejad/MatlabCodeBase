function Timelapse(row,col,site)
%row=4;col=5;site=5; %3_5_1 5_6_1
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\Documents\Projects\';
imagepath='E:\';
shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20140808 NLSpipCy2\20141024 GemVen mKatePIPcomp2\';

%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
datadir=[projectpath,experimentpath,'Data_Test3\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,shot,'\',shot,'_'];
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
SF=1;EF=157; %1:157
nucr=12;
maxjump=nucr*2; %default nucr*3
debrisarea=100;
boulderarea=1500; %MCF-10A: H2B:20 NLS:10
blobthreshold=-0.02;  %default 10xbin1=-0.02 20xbin2=-0.03; NLS=-0.01;
blurradius=3;
jitterwindow=3*nucr; %default 10xbin1=48
ringcalc=0; %set to zero if ring is unneeded
name1='H2B_'; %nuc
name2='Geminin_';
name3='mKate2_';
namenucedge='nucedge_';
adjustframe=0; %default
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
tracedata=ones(maxcellnum,EF,7+ringcalc*2)*NaN; %default 7
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
timetotal=tic;
%[firstgoodindex,blurthreshhigh,numthresh,badframes,height,width]=timelapsesetup(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_3(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
%[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_4thresh(rawdir,name1,frames,nucr,debrisarea,badframes,maskwrite);
%blurthreshlow=0;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([biasdir,name1,num2str(site),'.mat']); bias1=bias;
load([biasdir,name2,num2str(site),'.mat']); bias2=bias;
load([biasdir,name3,num2str(site),'.mat']); bias3=bias;
regheight=1:0.5*height; regwidth=1:0.5*width;
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw1=single(imread([rawdir,name1,num2str(f),'.tif'])); raw1=(raw1-bgcmos)./bias1;
    raw2=single(imread([rawdir,name2,num2str(f),'.tif'])); raw2=(raw2-bgcmos)./bias2;
    raw3=single(imread([rawdir,name3,num2str(f),'.tif'])); raw3=(raw3-bgcmos)./bias3;
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
    if i==firstgoodindex
        nuc_mask=threshmask(raw1,blurradius);
        nuc_mask=markershed(nuc_mask,round(nucr*2/3));
        foreground=nuc_mask;
        nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
        nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
        nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
    else
        nuc_mask=threshmask(raw1,1);
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
        jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        markerx=round(nuc_center(:,1)-jitters(f,1));
        markery=round(nuc_center(:,2)-jitters(f,2));
        markernan=isnan(markerx) | markerx<1 | markerx>width | markery<1 | markery>height;
        markerx(markernan)=[]; markery(markernan)=[];
        markeridx=sub2ind([height,width],markery,markerx);
        marker_mask=zeros(height,width);
        marker_mask(markeridx)=1; marker_mask=imdilate(marker_mask,strel('disk',4,0));
        [marker_label,markernum]=bwlabel(marker_mask);
        maskedmarkerlabels=marker_label.*nuc_mask;
        uniquemarkerlabels=unique(maskedmarkerlabels); uniquemarkerlabels(1)=[];
        missinglabels=find(~ismember(1:markernum,uniquemarkerlabels));
        marker_mask(ismember(marker_label,missinglabels))=0;
        nuc_mask=markershed_apriori(nuc_mask,marker_mask);
        foreground=nuc_mask;
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=imclearborder(nuc_mask);
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression=4;
    nanmask=imdilate(foreground,strel('disk',nucr/2));
    %nanmaskcyto=imdilate(foreground,strel('disk',nucr));
    %blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
    blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
    blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
    real1=bgsubmasked_global_2(raw1,nanmask,11,compression,50);
    real2=bgsubmasked_global_2(blur2,nanmask,11,compression,50);
    real3=bgsubmasked_global_2(blur3,nanmask,11,compression,50);
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)   
        fprintf('badframe: frame %0.0f\n',f);
        badframes(f)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(badextractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
            %imwrite(uint16(real2),[maskdir,'\',name2,num2str(f),'.tif']);
            %imwrite(uint16(real3),[maskdir,'\',name3,num2str(f),'.tif']);
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
        nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        %%% track & correct merges (update centers, masses and labels) %%%%
        H2BvsNLS=1; %1:H2B 2:NLS
        debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
        %[tracedata,curdata,tracking,nuc_label]=adaptivetrack_5(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),H2BvsNLS,debugpackage);
        %[tracedata,curdata,tracking,nuc_label]=adaptivetrack_6(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),H2BvsNLS,debugpackage);
        %[tracedata,curdata,tracking,nuc_label]=adaptivetrack_7_winrad(lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),maxjump,H2BvsNLS,debugpackage);
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_8_split(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),maxjump,debrisarea,H2BvsNLS,debugpackage);
        badframes(f)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
        %imwrite(uint16(real2),[maskdir,'\',name2,num2str(f),'.tif']);
        %imwrite(uint16(real3),[maskdir,'\',name3,num2str(f),'.tif']);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=mean(real1(nuc_info(cc).PixelIdxList));
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
        sig3(n)=median(real3(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        innerrad=1; outerrad=5; %10xB1|20xB2: 1/5
        ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
        ring_info=regionprops(ring_label,'PixelIdxList');
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_fgmode=nanvec;
        for n=1:numlivecells
            cc=cellid(n);
            ring2all=real2(ring_info(cc).PixelIdxList);
            sig2ring_75th(n)=prctile(ring2all,75);
            ring2foreground=ring2all(ring2all>0);
            if numel(ring2foreground)<50
                 ring2foreground=ring2all;
            end
            if numel(ring2all)<50
                 continue;
            end             
            sig2ring_fgmedian(n)=nanmedian(ring2foreground);
            numbins=25;
            sig2ring_fgmode(n)=getmode(ring2foreground,numbins);
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc==1
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode];
    else
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3];
    end
    if maxcellnum-max(cellid)<blocksize
        tempdata=ones(blocksize,EF,7+ringcalc*2)*NaN; %default 7
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
%[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr,H2BvsNLS);
[tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
toc(timetotal);
clear all; clear mex;

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=imadjust(mat2gray(raw1));
tempframe(:,:,3)=marker_mask;
tempframe(:,:,1)=extractmask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}