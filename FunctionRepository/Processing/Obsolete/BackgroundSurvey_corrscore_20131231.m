function Immunostain_2(rowmat,colmat,sitemat)
rowmat=1:7;colmat=1:12;sitemat=1:4;
%rowmat=1;colmat=3;sitemat=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath='2013-12-12_BackgroundCharacterization\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
numsites=numel(rowmat)*numel(colmat)*numel(sitemat);
%medint=ones(numsites,1)*NaN;
meanint=ones(numel(rowmat),numel(colmat),numel(sitemat))*NaN;
stdint=meanint; relint=meanint;
cc=0;
time1=tic;
for row=rowmat
    for col=colmat
        for site=sitemat
            cc=cc+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
            %[medint(cc),stdint(cc),relint(cc)]=main(rawdir);
            [meanint(row,col,site),stdint(row,col,site),relint(row,col,site)]=main(rawdir);
        end
    end
end
scatter(meanint(:),stdint(:));
axis([600 700 0 50]);
emedint=meanint; %emedint(22)=[];
estdint=stdint; %estdint(22)=[];
erelint=relint; %erelint(22)=[];
toc(time1);

fun=@(block_struct) nanmean(block_struct.data(:));  %function handle for blockproc call
a=single(imread('H:\Images\2013-12-12_BackgroundCharacterization\Raw\1_1_1_Cy3_stain.tif'));
%a=single(imread('H:\Images\2013-11-26_CycDp21SerumRelease\Experiment_20131216\Raw\1_3_1_Cy3_stain.tif'));
[height,width]=size(a);
blocknum=10;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
a=blockproc(a,[blockheight blockwidth],fun);
%a=imresize(a,0.1);
b=single(imread('H:\Images\2013-12-12_BackgroundCharacterization\Raw\7_12_2_Cy3_stain.tif'));
%b=single(imread('H:\Images\2013-11-26_CycDp21SerumRelease\Experiment_20131216\Raw\1_11_1_Cy3_stain.tif'));
b=blockproc(b,[blockheight blockwidth],fun);
%a=a(:); meana=mean(a); stda=std(a); norma=(a-meana)/stda; b=b(:); meanb=mean(b); stdb=std(b); normb=(b-meanb)/stdb; axb=norma.*normb; c=sum(axb(:))/(numel(a));
a=a(:); a=a-mean(a); b=b(:); b=b-mean(b); c=dot(a/norm(a),b/norm(b));
keyboard;
%imagesc(c); colorbar;
%set(gcf,'color','w','PaperPosition',[0 0 4 3]);
%saveas(gcf,'h:\Downloads\Fig.jpg');
end



function [meanint,stdint,relint]=main(rawdir)
name='Cy3_';
raw=imfilter(single(imread([rawdir,name,'stain.tif'])),fspecial('disk',10));
rmask=raw>900; rmask=imdilate(rmask,strel('disk',100));
raw(rmask)=NaN;
%raw(raw>1000)=NaN;

% rawsmall=imresize(raw,0.01);
% 
% medint=median(raw(:));
% stdint=std(raw(:));
% relint=abs(raw(:)-medint)/medint;
% relint=median(relint);

fun=@(block_struct) nanmedian(block_struct.data(:));  %function handle for blockproc call
[height,width]=size(raw);
blocknum=10;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
rawblock=blockproc(raw,[blockheight blockwidth],fun);
meanint=mean(rawblock(:));
stdint=std(rawblock(:));
relint=abs(rawblock(:)-meanint)/meanint;
relint=median(relint);
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