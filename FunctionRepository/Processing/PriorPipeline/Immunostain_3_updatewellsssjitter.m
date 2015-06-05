function Immunostain_3_updatewellsssjitter(row,col,site)
%row='D';col='05';site='4';
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\';
datadir = ([path,'Data\']);
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
rawdir = [path,'Raw\',shot,'\'];
cpdir = [path,'Processed\',shot];

nucr = 12;
nucname = 'CFP_'; %nuc
nucedgename = 'nucedge&ring_';
YFPname = 'YFP_'; %DHB
RFPname = 'TexasRed_'; %p21
IFname = 'Cy5_';
DAPIname = 'DAPI_';

debrisarea=pi*(nucr/4)^2;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucfile=[rawdir,nucname,'poststain.tif'];
nuc_raw=single(imread(nucfile));
[height,width]=size(nuc_raw);
blocknum=3;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);

timetotal=tic;
%%% copy images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NEfile=[cpdir,'\',nucedgename,'poststain.tif'];
YFPfile=[rawdir,YFPname,'poststain.tif'];
RFPfile=[rawdir,RFPname,'poststain.tif'];
IFfile=[rawdir,IFname,'poststain.tif'];
DAPIfile=[rawdir,DAPIname,'poststain.tif'];

%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YFP_raw=single(imread(YFPfile));
RFP_raw=single(imread(RFPfile));
IF_raw_post=single(imread(IFfile));
DAPI_raw=single(imread(DAPIfile));

%%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_bs=bgsub(log(nuc_raw),10*nucr,0.05);
YFP_bs=bgsub_MC(log(YFP_raw),blockheight,blockwidth);
RFP_bs=bgsub(log(RFP_raw),10*nucr,0.05);
IF_bs_post=bgsub(log(IF_raw_post),10*nucr,0.05);
DAPI_bs=bgsub(log(DAPI_raw),10*nucr,0.05);

%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calcbleedthroughrate = 1; %default method
prestaincorrection = 0;
if calcbleedthroughrate==1
    load([datadir,'bleedthroughrate.mat'],'bleedthroughrate');
    IF_bleedthrough=RFP_bs*bleedthroughrate(1)+bleedthroughrate(2);
    IF_bs=IF_bs_post-IF_bleedthrough;
elseif prestaincorrection==1
    nuc_raw_pre=single(imread([rawdir,nucname,'prestain.tif']));
    IF_raw_pre=single(imread([rawdir,IFname,'prestain.tif']));
    IF_bs_pre=bgsub(log(IF_raw_pre),10*nucr,0.05);
    [ppjx,ppjy]=CalcJitter(nuc_raw_pre,nuc_raw);
    ppjx=round(ppjx); ppjy=round(ppjy);
    prey=1+ppjy*(ppjy>0):height+ppjy*(ppjy<0);
    prex=1+ppjx*(ppjx>0):width+ppjx*(ppjx<0);
    posty=1-ppjy*(ppjy<0):height-ppjy*(ppjy>0);
    postx=1-ppjx*(ppjx<0):width-ppjx*(ppjx>0);
    IF_bs_pre_shifted=zeros(height,width);
    IF_bs_pre_shifted(posty,postx)=IF_bs_pre(prey,prex);
    IF_bs=IF_bs_post-IF_bs_pre_shifted;
    nuc_bs_temp=nuc_bs;
    nuc_bs=zeros(height,width);
    nuc_bs(posty,postx)=nuc_bs_temp(posty,postx); %remove nuclei from crop margins
else
    IF_bs=IF_bs_post;
end

%%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
nuc_mask=getnucmask_histanddeflection(nuc_bs,nucr,1,1);
[nuc_info,nuc_label,finalcytoring]=simplercytoring(nuc_mask,nucr);
ring_info=regionprops(finalcytoring,'PixelIdxList');
%%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(nuc_info,1);
numrings=size(ring_info,1);
if numcells>numrings
    nucexclude=[numrings+1:numcells];
    nuc_info(nucexclude)=[];
end
nucexclude=zeros(numrings,1);
for k=1:numrings
    if nuc_info(k).Area < debrisarea
        nucexclude(k)=1;
    end
end
nucexclude=find(nucexclude);
nuc_info(nucexclude)=[];
ring_info(nucexclude)=[];
numcells=size(nuc_info,1);

%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label(ismember(nuc_label,nucexclude))=0;
finalcytoring(ismember(finalcytoring,nucexclude))=0;
%extractmask=bwmorph(nuc_label,'remove') + logical(finalcytoring);
extractmask=bwmorph(nuc_label,'remove');
imwrite(uint16(extractmask),NEfile);
%%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tv=zeros(numcells,1);
XX=tv; YY=tv; nuc_area=tv;
nuc_median=tv;
YFP_nucmedian=tv;
RFP_median=tv;
YFP_ringmean=tv;
IF_median=tv;
DAPI_mean=tv;
for cc=1:numcells
    XX(cc,1)=nuc_info(cc).Centroid(1);  %x value of centroid
    YY(cc,1)=nuc_info(cc).Centroid(2);  %y value of centroid
    nuc_area(cc,1)=nuc_info(cc).Area;

    nuc_median(cc,1)=median(nuc_bs(nuc_info(cc).PixelIdxList));
    RFP_median(cc,1)=median(RFP_bs(nuc_info(cc).PixelIdxList));

    YFP_nucmedian(cc,1)=median(YFP_bs(nuc_info(cc).PixelIdxList));
    %YFP_ringmedian(cc,1)=median(YFP_bs(ring_info(cc).PixelIdxList));
    allringpixels=YFP_bs(ring_info(cc).PixelIdxList);
    topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
    YFP_ringmean(cc,1)=mean(topringpixels);
    
    IF_median(cc,1)=median(IF_bs(nuc_info(cc).PixelIdxList));
    DAPI_mean(cc,1)=mean(DAPI_bs(nuc_info(cc).PixelIdxList));
end
IFdata=[XX,YY,nuc_median,nuc_area,YFP_nucmedian,RFP_median,YFP_ringmean,IF_median,DAPI_mean]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame

%%% update wellsss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'wellsss_',shot,'.mat'],'wellsss','badframes');
lastframe = size(wellsss,3); %last frame of normal timelapse
stainframe = lastframe+1;
for f=1:lastframe
    temparray=wellsss{f};
    numrows=size(temparray,1);
    wellsss{f}=cat(2,wellsss{f},zeros(numrows,2));
end
wellsss{1,1,stainframe}=IFdata; %appends immunostain frame
save([datadir,'wellsss_IF_',shot,'.mat'],'wellsss','badframes');

%%% calculate last jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lastframe_mask=single(imread([cpdir,'\',nucedgename,num2str(lastframe),'.tif']));
lastframe_mask=imfill(lastframe_mask);
stain_mask=single(extractmask);
stain_mask=imfill(stain_mask);
[lastjitterx,lastjittery]=get_alignment_shift(lastframe_mask,stain_mask,100);

%%% update jitter list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'jitter_',shot,'.mat'],'x','y');
%%% one-time correction for 20130715 set %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lastjitterx=x(end)-lastjitterx;   %subtract for get_alignment_shift
lastjittery=y(end)-lastjittery;
x=[x lastjitterx];
y=[y lastjittery];
save([datadir,'jitter_IF_',shot,'.mat'],'x','y');

toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(REs_bs));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}
cd([codepath,'Processing\']); %return to this directory