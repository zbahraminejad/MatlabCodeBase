function X_immunofluorescence_simple(row)
%cd('h:\Documents\Timelapse\Code\Development\Functions'); %change directory for function calls
path = 'h:\Documents\Immunofluorescence\2012-12-11 PD0325901 Dose Response\';

if ismember(row,[2 4 6])
    experiment = 'pFra-1-488 pErk-594\';
else
    experiment = 'pErk-488 pAkt-594\';
end
rawdir = [path,experiment,'Raw\'];

cols = [1:12];
sites = [1:9];
%cols = [1];
%sites = [9];

%numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
totalcellcount=zeros(numcols,numsites);
totalmedianred=zeros(numcols,numsites);
totalmedianfitc=zeros(numcols,numsites);

nucroot = '_DAPI_';
redroot = '_Texas Red_';
fitcroot = '_FITC_';
nucr=12;
frame=1;


for colidx=1:numcols
    for siteidx=1:numsites
col=cols(colidx); site=sites(siteidx);
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
fprintf('row %0.0f, col %0.0f, site %0.0f\n',row,col,site);

%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_or=single(imread([rawdir,shot,nucroot,num2str(frame),'.tif']));
REs_or=single(imread([rawdir,shot,redroot,num2str(frame),'.tif']));
CEs_or=single(imread([rawdir,shot,fitcroot,num2str(frame),'.tif'])); 

%%% image processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_bl=log(imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric')); %take log to decrease the variance of the signal.  use fspecial to create a disk-shaped filter to make a blurred (bl) image
DAs_bs=bgsub(DAs_bl,10*nucr,0.05);  %background subtraction on the blurred image. search 10x the radius, sort, find the value of the bottom 5th percentile. can change these #s
REs_bl=(imfilter(REs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
REs_bs=bgsub(REs_bl,10*nucr,0.05);
CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
CEs_bs=bgsub(CEs_bl,10*nucr,0.05);

%%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_ma=getdapimask(DAs_bs,nucr);
[realnuc_la,numcells]=bwlabel(DAs_ma);
meandapi=cell2mat(struct2cell(regionprops(realnuc_la,DAs_bs,'MeanIntensity'))');
meanred=cell2mat(struct2cell(regionprops(realnuc_la,REs_bs,'MeanIntensity'))');
meanfitc=cell2mat(struct2cell(regionprops(realnuc_la,CEs_bs,'MeanIntensity'))');

totalcellcount(col,site)=numcells;
totalmediandapi(col,site)=median(meandapi);
totalmedianred(col,site)=median(meanred);
totalmedianfitc(col,site)=median(meanfitc);
    end
end
save([path,experiment,'R',num2str(row),'_alldata.mat'],'totalcellcount','totalmediandapi','totalmedianred','totalmedianfitc');
%cd('h:\Documents\Immunofluorescence\Code\');    %return to original directory
end