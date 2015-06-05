function Immunostain_2(rowmat,colmat,sitemat)
%rowmat=1:7;colmat=1:12;sitemat=1;
rowmat=1;colmat=3;sitemat=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath='2013-12-12_BackgroundCharacterization\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
numsites=numel(rowmat)*numel(colmat)*numel(sitemat);
medint=ones(numsites,1)*NaN; stdint=medint; relint=medint;
cc=0;
time1=tic;
for row=rowmat
    for col=colmat
        for site=sitemat
            cc=cc+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
            [medint(cc),stdint(cc),relint(cc)]=main(rawdir);
        end
    end
end
toc(time1);

fun=@(block_struct) nanmean(block_struct.data(:));  %function handle for blockproc call
a=single(imread('H:\Images\2013-12-12_BackgroundCharacterization\Raw\4_4_1_Cy3_stain.tif'));
%a=single(imread('H:\Images\2013-11-26_CycDp21SerumRelease\Experiment_20131216\\Raw\1_3_1_Cy3_stain.tif'));
[height,width]=size(a);
blocknum=10;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
a=blockproc(a,[blockheight blockwidth],fun);
%a=imresize(a,0.1);
b=single(imread('H:\Images\2013-12-12_BackgroundCharacterization\Raw\4_11_1_Cy3_stain.tif'));
%b=single(imread('H:\Images\2013-11-26_CycDp21SerumRelease\Experiment_20131216\\Raw\1_11_1_Cy3_stain.tif'));
b=blockproc(b,[blockheight blockwidth],fun);
%b=imresize(b,0.1);
c=a-b; c=c-median(c(:));
%cliplow=-100; cliphigh=100;
cliplow=645; cliphigh=695;
%cliplow=280; cliphigh=340;
a(a>cliphigh)=cliphigh; a(a<cliplow)=cliplow;
b(b>cliphigh)=cliphigh; b(b<cliplow)=cliplow;
c(c>25)=25;
imagesc(c); colorbar;
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end



function [medint,stdint,relint]=main(rawdir)
name='Cy3_';
raw=single(imread([rawdir,name,'stain.tif']));

fun=@(block_struct) nanmedian(block_struct.data(:));  %function handle for blockproc call
[height,width]=size(raw);
blocknum=20;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
rawblock=blockproc(raw,[blockheight blockwidth],fun);
medint=median(rawblock(:));
stdint=std(rawblock(:));
relint=abs(rawblock(:)-medint)/medint;
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