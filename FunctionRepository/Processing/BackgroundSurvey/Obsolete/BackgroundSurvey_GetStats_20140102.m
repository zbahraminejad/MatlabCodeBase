function BackgroundSurvey_getstats
rowmat=1:7;colmat=1:12;sitemat=1:4;
%rowmat=1;colmat=3;sitemat=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
experimentpath='2013-12-12_BackgroundCharacterization\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
filetag='Cy3_';
meanint=ones(numel(rowmat),numel(colmat),numel(sitemat))*NaN; stdint=meanint;
cc=0;
for row=rowmat
    for col=colmat
        for site=sitemat
            cc=cc+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            rawdir = [imagepath,experimentpath,'Raw\'];
            [meanint(row,col,site),stdint(row,col,site)]=getstats(rawdir,shot,filetag);
        end
    end
end
scatter(meanint(:),stdint(:));
%imagesc(1:12,1:7,meanint(:,:,1)); colorbar;
xlabel('mean RFU'); ylabel('stddev RFU');
set(gcf,'color','w','PaperPosition',[0 0 8 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end

function [meanint,stdint]=getstats(rawdir,shot,filetag)
filename=[rawdir,shot,'_',filetag,'stain.tif'];
raw=imfilter(single(imread(filename)),fspecial('disk',10));
rmask=raw>900; rmask=imdilate(rmask,strel('disk',100));
raw(rmask)=NaN;
fun=@(block_struct) nanmedian(block_struct.data(:));  %function handle for blockproc call
[height,width]=size(raw);
blocknum=10;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
rawblock=blockproc(raw,[blockheight blockwidth],fun);
meanint=mean(rawblock(:));
stdint=std(rawblock(:));
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}