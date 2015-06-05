%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagepath='D:\Images\';
imagepath='G:\Zahra\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='150507_LiveCell\';
% staindir = [experimentpath,'Stains\'];
 biasdir=[imagepath,experimentpath,'Raw\Bias\'];
%biasdir=[imagepath,experimentpath,'Raw\','Bias\'];
separatedirectories=0;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucradius=8;
names={
%     'hoechst';
    'Cy5';
};
rowmat=2:7;
colmat=2:11;
sitemat=1:4;
framemat=[1 100 200 300 400 500];
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat); numframes=numel(framemat);
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos=BG; 
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos=BG; 
% load([shadingpath,'cameranoise_4sites.mat'],'smoothnoise'); bgcmos=smoothnoise; % 1836x1836
% bgcmos = imresize(bgcmos,0.5); % for bin 2 images 1080x1080
[height,width]=size(bgcmos);
nucr=8;
debrisarea=50;
blobthreshold=-0.03;
for i=1:length(names)
    for s=1:numsites
        site=sitemat(s);
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
                        rawdir=[imagepath,experimentpath,'Raw\Day4\'];
                        
                        %raw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        %raw=single(imread([rawdir,shot,'\',names{i},'_',num2str(frame),'.tif']));
                        nucraw=double(imread([rawdir,shot,'_',names{i},'.tif']));
                    else
                        % use this for regular time lapse stuff
                        rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
                        nucraw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        % use this for fixed cell measurements of bias
%                         rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
%                         rawdir=[imagepath,staindir,shot,'_'];
%                         nucraw=double(imread([rawdir,names{i},'.tif']));
                        
                    end
                    nucraw=nucraw-bgcmos;
%                     nucmask=blobdetector_3_bin2(log(nucraw),nucr,blobthreshold,debrisarea);
                    nucmask=threshmask(nucraw,3);
                    nucmask=imdilate(nucmask,strel('disk',round(nucradius/2),0));

                    blocknum=31;
                    prctilethresh=10;
                    blurnan=imfilter(nucraw,fspecial('disk',3),'symmetric'); %10x:3, 10xbin2:2, 20x:6
                    blurnan(nucmask)=NaN;
                    %bg=blockpercentile(blurnan,blocknum,prctilethresh);
                    bgblock=blockpercentile_blockimage(blurnan,blocknum,prctilethresh);

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
        blockbias=nanmedian(biasstack,3);
        bias=imresize(blockbias,[height width],'bicubic');
        %%%%for immunostain
%         save([biasdir,names{i},'_Stain_',num2str(site),'.mat'],'bias');
%         save([biasdir,names{i},'_',num2str(site),num2str(1),'.mat'],'bias');
        %%%%for regular timelapse bias
        save([biasdir,names{i},'_',num2str(site),'.mat'],'bias');
    end
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nucmask,'remove');
[height,width]=size(nucmask);
RGB=zeros(height,width,3);
RGB(:,:,1)=imadjust(mat2gray(nucraw));
RGB(:,:,2)=extractmask;
figure,imshow(RGB);
%}