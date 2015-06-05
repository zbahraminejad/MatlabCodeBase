%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagepath='D:\Images\';
%imagepath='D:\Images\';
imagepath='G:\Zahra\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='150507_LiveCell\';
staindir = [experimentpath,'Stains\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=0;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucradius=8;
names={
%     'mCherry';
'YFP';

%       'CFP';
%       'cebpb';
};
nucname = 'Cy5';
rowmat=2:7;
colmat=2:11;
sitemat=1:4;
framemat=[ 1 100 200 300 400 500];%
maskdilation=round(nucradius*round(1.5));
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat); numframes=numel(framemat);
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG')r; bgcmos=BG; 
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG;
% bgcmos = imresize(bgcmos,0.5); % for bin 2 images 1080x1080
% load([shadingpath,'cameranoise_4sites.mat'],'smoothnoise'); bgcmos=smoothnoise; % 1836x1836
[height,width]=size(bgcmos);
nucr=8;
debrisarea=25;
blobthreshold=-0.03;
for i=1:length(names)
    for s=1:numsites
        site=sitemat(s);
        %load([Biasdir,nucname,'_',num2str(site),'.mat'],'bias');
%         load([biasdir,nucname,'_Stain_',num2str(site),'.mat'],'bias');
        load([biasdir,nucname,'_',num2str(site),'.mat'],'bias');
        nucbias=bias;
        biasstack=[];
        for c=1:numcols
            col=colmat(c);
            for r=1:numrows
                row=rowmat(r);
                for f=1:numframes
                    frame=framemat(f);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    wellname = nameandsite(shot);
                    if separatedirectories
                        rawdir=[imagepath,experimentpath,'Raw\'];
                        %raw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        %raw=single(imread([rawdir,shot,'\',names{i},'_',num2str(frame),'.tif']));
                        nucraw=single(imread([rawdir,shot,'_',nucname,'.tif']));
                        sigraw=single(imread([rawdir,shot,'_',names{i},'.tif']));
                    else
                        %use this for regular time lapse calculuations
                        rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
                        nucraw=double(imread([rawdir,nucname,'_',num2str(frame),'.tif']));
                        sigraw=double(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        %use this for fixed cell imaging
%                         rawdir=[imagepath,staindir,shot,'_'];
%                         nucraw=single(imread([rawdir,nucname,'.tif']));
%                         sigraw=single(imread([rawdir,names{i},'.tif']));
                    end
                    nucraw=(nucraw-bgcmos)./nucbias;
%                     nucmask=threshmask(nucraw,3);
                    nucmask=blobdetector_3_bin2(log(nucraw),nucr,blobthreshold,debrisarea);
                    
                    % Blur signal image
                    sigblur=imfilter(sigraw,fspecial('gaussian',round(nucradius/2)),'symmetric');
                    % Subtract camera noise from image
                    backgroundonly=sigblur-bgcmos;
                    % Remove foreground pixels by setting them to NaN (Not a Number)
                    nanmask=imdilate(nucmask,strel('disk',maskdilation,0));
                    backgroundonly(nanmask)=NaN;
                    
                    % Break the image into sections and estimate the background in each section
                    % by a given intensity percentile. Interpolate all of the sections to
                    % return a smoothened background image.
                    blocknum=31;
                    prctilethresh=10;
                    %background=blockpercentile(backgroundonly,blocknum,prctilethresh);
                    bgblock=blockpercentile_blockimage(backgroundonly,blocknum,prctilethresh);
                    
                    midrc=ceil(blocknum/2);
                    refval=bgblock(midrc);
                    if ~isnan(refval)
                        %bgblocknorm=bgblock/max(bgblock(:));
                        bgblocknorm=bgblock/refval;
                        biasstack=cat(3,biasstack,bgblocknorm);
                    end
                end
            end
        end
        %bias=median(biasstack,3);
        blockbias=nanmedian(biasstack,3);
        bias=imresize(blockbias,[height width],'bicubic');
%         save([biasdir,names{i},'_',num2str(site),'.mat'],'bias');
        save([biasdir,names{i},'_',num2str(site),'.mat'],'bias');
%         save([biasdir,names{i},'_Stain_',num2str(site),'.mat'],'bias');
    end
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nanmask,'remove');
[height,width]=size(nucmask);
RGB=zeros(height,width,3);
RGB(:,:,1)=imadjust(mat2gray(sigraw));
RGB(:,:,2)=extractmask;
figure,imshow(RGB);
%}