function Timelapse_bin2(row,col,site)
% row=7;col=6;site=3;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projectpath='D:\Documents\Projects\';
% imagepath='G:\Projects\CellCycle_Differentiation\';
% imagepath = 'E:\Michael\';
imagepath = 'G:\Zahra\';
%imagepath='E:\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';swe
% experimentpath='12052014-CellCycle-24hr-48hr-bins\';
% experimentpath= '20150126-gem-pparg-diff\';
experimentpath = '150507_LiveCell\';
%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
wellname = nameandsite(shot);
datadir=[imagepath,experimentpath,'Data\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    
    %rawdir=[imagepath,experimentpath,shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot];
    %maskdir=[imagepath,experimentpath,'Mask\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\'];
%     realdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
maskwrite=1;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
if separatedirectories==0
    maskdir=[maskdir,shot,'_'];
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SF=1;EF=584; %1:319
SF = 1;EF = 555;
nucr=8;
maxjump=nucr*3; %default nucr*3
debrisarea=25;
boulderarea=1000; %MCF-10A: H2B:20 NLS:10
blobthreshold=-0.03;  %default 10xbin1=-0.02 20xbin2=-0.03; NLS=-0.01;
%blurradius=3;
ringcalc=0; %set to zero if ring is unneeded
name1='Cy5_'; %nuc
name2='YFP_';
% name3='YFP_';
% name4='mCherry_';
namenucedge='nucedge_';
% adjustframe=[114,386]; %default
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
framelength = length(SF:EF);
totalframes=numel(frames);
badframes=ones(framelength,1)*NaN;
if SF>1
    badframes(1:SF-1)=0;
end
jitters=zeros(EF,2);
blocksize=10000;
maxcellnum=blocksize;
if ringcalc==1
    parameternum=11; % 13 for two rings
    tracedata=ones(maxcellnum,EF,parameternum)*NaN;
else
    parameternum=6; % 6 for 1 signal 7 for 2 signals
    tracedata=ones(maxcellnum,EF,parameternum)*NaN;
end
tracking=ones(maxcellnum,5)*NaN; %[mother,mergingcell1,mergingcell2,mergestart,mergeend]
cellcounter = ones(framelength,1);
timetotal=tic;
%[firstgoodindex,blurthreshhigh,numthresh,badframes,height,width]=timelapsesetup(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_3(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
% [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_4thresh(rawdir,name1,frames,nucr,debrisarea,badframes,maskwrite);
%blurthreshlow=0;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos=BG;
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG; 
load([biasdir,name1,num2str(site),'.mat']); bias1=bias;
load([biasdir,name2,num2str(site),'.mat']); bias2=bias;
% load([biasdir,name3,num2str(site),'.mat']); bias3=bias;
% load([biasdir,name4,num2str(site),'.mat']); bias4=bias;
regheight=1:0.5*height; regwidth=1:0.5*width;
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    
    timeframe=tic;
        
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    raw1=single(imread([rawdir,name1,num2str(f),'.tif'])); raw1=(raw1-bgcmos)./bias1;
    raw2=single(imread([rawdir,name2,num2str(f),'.tif'])); raw2=(raw2-bgcmos)./bias2;
%     raw3=single(imread([rawdir,name3,num2str(f),'.tif'])); raw3=(raw3-bgcmos)./bias3;
%     raw4=single(imread([rawdir,name4,num2str(f),'.tif'])); raw4=(raw4-bgcmos)./bias4;
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==firstgoodindex
        %%% LoG filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=blobdetector_3_bin2(log(raw1),nucr,blobthreshold,debrisarea);
        foreground=nuc_mask;
%         %%% Double-Threshold filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=threshmask(raw1,blurradius);
%         nuc_mask=markershed(nuc_mask,round(nucr*2/3));
%         foreground=nuc_mask;
%         nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=bwareaopen(nuc_mask,debrisarea);
        %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
        nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
    else
%         nuc_mask=threshmask(raw1,1);
        nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
        jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        %%% LoG filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
%         %%% Double-Threshold filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         nuc_mask=threshmask(raw1,blurradius);
%         nuc_mask=markershed(nuc_mask,round(nucr*2/3));
%         foreground=nuc_mask;
%         nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea*2);
%         %%% a priori marker-based watershed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [nuc_mask,marker_mask]=aprior_markermask(nuc_mask,nuc_center,jitters(f,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        foreground=nuc_mask;
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask=imclearborder(nuc_mask);
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression1=2; % compression 4 for 2160x2160 and 2 for 1080x1080
%     compression2=1;
    nanmask=imdilate(foreground,strel('disk',round(nucr*1.5)));
    %nanmaskcyto=imdilate(foreground,strel('disk',nucr));
    %blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
    blur2=imfilter(raw2,fspecial('disk',3),'symmetric'); % default is blur by 3 for 10x bin1, but 2 for 10x bin2
%     blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
%     blur4=imfilter(raw4,fspecial('disk',3),'symmetric');

    real1=bgsubmasked_global_2(raw1,nanmask,1,compression1,50);    
    real2=bgsubmasked_global_2(blur2,nanmask,1,compression1,50);
%     real3=bgsubmasked_global_2(blur3,nanmask,1,compression2,50);
%     real4=bgsubmasked_global_2(blur4,nanmask,1,compression1,50);
    
%     if (row==4|| row==5)&&(col==2||col==3) 
%     imwrite(uint16(real1),[realdir,name1,num2str(f),'.tif']);
%     imwrite(uint16(real2),[realdir,name2,num2str(f),'.tif']);
%     imwrite(uint16(real3),[realdir,name3,num2str(f),'.tif']);
%     imwrite(uint16(real4),[realdir,name4,num2str(f),'.tif']);
%     end
    
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
            imwrite(uint16(badextractmask),[maskdir,namenucedge,num2str(f),'.tif']);
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
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_8_split(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,nucr,jitters(f,:),maxjump,debrisarea*2,H2BvsNLS,debugpackage);
        badframes(f)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    %nuc_inner=imerode(nuc_label,strel('disk',1,0));
    %extractmask=nuc_label-nuc_inner;
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,namenucedge,num2str(f),'.tif']);
        %imwrite(uint16(real2),[maskdir,name2,num2str(f),'.tif']);
        %imwrite(uint16(real3),[maskdir,name3,num2str(f),'.tif']);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    cellcounter(f) = numlivecells;
    disp(numlivecells)
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; %sig3=nanvec; sig4=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=mean(real1(nuc_info(cc).PixelIdxList));
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
%         sig3(n)=median(real3(nuc_info(cc).PixelIdxList));
%         sig4(n)=median(real4(nuc_info(cc).PixelIdxList));
    end
    if ringcalc==1
        innerrad=1; outerrad=5; %10xB1|20xB2: 1/5
        ring_label=getcytoring_thicken(nuc_label,innerrad,outerrad,real2);
        ring_info=regionprops(ring_label,'PixelIdxList');
%         bgthreshold = abs(prctile(real2(~nanmask),1)); % use only with median bg subtraction, not with mode bg subtraction
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_fgmode=nanvec;
%         sig3ring_75th=nanvec; sig3ring_fgmedian=nanvec; sig3ring_fgmode=nanvec;
        for n=1:numlivecells
            cc=cellid(n);
            if cc>numel(ring_info)
                break;
            end
            ring2all=real2(ring_info(cc).PixelIdxList);
             ring2all(ring2all>prctile(ring2all,98))=[];
%             figure,histogram(ring2all)
            sig2ring_75th(n)=prctile(ring2all,75);
            ring2foreground=ring2all(ring2all>0);
            if numel(ring2foreground)<30 %% was 100
                 ring2foreground=ring2all;
            end
%             ring3all=real3(ring_info(cc).PixelIdxList);
%             ring3all(ring3all>prctile(ring3all,98))=[];
%             sig3ring_75th(n)=prctile(ring3all,75);
%             ring3foreground=ring3all(ring3all>0);
%             if numel(ring3foreground)<100
%                  ring3foreground=ring3all;
%             end
            if numel(ring2all)>50
                 sig2ring_fgmedian(n)=nanmedian(ring2foreground);
                 numbins=15;
                 sig2ring_fgmode(n)=getmode(ring2foreground,numbins);
%                  sig3ring_fgmedian(n)=nanmedian(ring3foreground);
%                  sig3ring_fgmode(n)=getmode(ring3foreground,numbins);
            end
        end
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ringcalc == 1
%         tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode,sig3,sig3ring_75th,sig3ring_fgmedian,sig3ring_fgmode];
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgmode,sig3,sig4];
    else 
%         tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3];  
        tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2]; 
    end
    if maxcellnum-max(cellid)<blocksize
        tempdata=ones(blocksize,EF,parameternum)*NaN;
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
% [tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%[tracedata,genealogy,jitters]=postprocessing(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr,H2BvsNLS);
[tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',shot,'_nolink','.mat'],'tracedata','genealogy','jitters');
toc(timetotal);
clear all; 
clear mex;

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}