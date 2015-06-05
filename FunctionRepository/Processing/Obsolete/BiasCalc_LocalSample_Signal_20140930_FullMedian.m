%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath='D:\Images\';
shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20130607 p21dCy2\20140208 H2B DHB p21dCy1dK\';
Biasdir=[imagepath,experimentpath,'Raw\Bias\'];
separatedirectories=1;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucname='H2B';
names={
    %'DHB';
    %'p21dCy1dK';
    'p21';
};
rowmat=[3:6]; %[3 4 7]
colmat=[7]; %[5:7]
sitemat=1:4; %1
framemat=[1];
nucradius=12;
maskdilation=round(nucradius/2);
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat); numframes=numel(framemat);
load([shadingpath,'BG.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
for i=1:length(names)
    for s=1:numsites
        site=sitemat(s);
        load([Biasdir,nucname,'_',num2str(site),'.mat'],'bias');
        nucbias=bias;
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
                        nucraw=single(imread([rawdir,shot,'\',nucname,'_stain.tif']));
                        sigraw=single(imread([rawdir,shot,'\',names{i},'_stain.tif']));
                    else
                        rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
                        nucraw=single(imread([rawdir,nucname,'_stain.tif']));
                        sigraw=single(imread([rawdir,names{i},'_stain.tif']));
                    end
                    nucraw=(nucraw-bgcmos)./nucbias;
                    nucblur=imfilter(nucraw,fspecial('gaussian',round(nucradius/2)),'symmetric');
                    normlograw1=mat2gray(log(nucblur));
                    T=graythresh(normlograw1);
                    nucmask=im2bw(normlograw1,T);
                    nucmask=imfill(nucmask,'holes');
                    
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
                    blocknum=5;
                    prctilethresh=10;
                    background=blockpercentile(backgroundonly,blocknum,prctilethresh);
                    maxillumination=max(background(:));
                    bgnorm=background/maxillumination;
                    
                    biasstack=cat(3,biasstack,bgnorm);
                end
            end
        end
        bias=median(biasstack,3);
        save([Biasdir,names{i},'_',num2str(site),'.mat'],'bias');
        %save([imagepath,experimentpath,'Raw\IBG_',names{i},'.mat'],'bias');
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