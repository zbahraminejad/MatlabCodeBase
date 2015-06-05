%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath='D:\Images\';
shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20130607 p21dCy2\20140208 H2B DHB p21dCy1dK\';
Biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=1;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucradius=12;
names={
    'H2B';
};
rowmat=[3:6]; %[3 4 7]
colmat=[7]; %[5:7]
sitemat=1:4; %1
framemat=[1];
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat); numframes=numel(framemat);
load([shadingpath,'BG.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
[height,width]=size(bgcmos);
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
                    if separatedirectories
                        rawdir=[imagepath,experimentpath,'Raw\'];
                        %raw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        %raw=single(imread([rawdir,shot,'\',names{i},'_',num2str(frame),'.tif']));
                        nucraw=single(imread([rawdir,shot,'\',names{i},'_stain.tif']));
                    else
                        rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
                        nucraw=single(imread([rawdir,names{i},'_stain.tif']));
                    end
                    nucraw=nucraw-bgcmos;
                    nucmask=threshmask(nucraw,3);
                    nucmask=imdilate(nucmask,strel('disk',round(nucradius/2),0));

                    blocknum=9;
                    prctilethresh=10;
                    blurnan=imfilter(nucraw,fspecial('disk',3),'symmetric'); %10x:3 20x:6
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
        save([Biasdir,names{i},'_',num2str(site),'.mat'],'bias');
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