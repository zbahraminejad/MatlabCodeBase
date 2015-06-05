function Immunostain_2(rowmat,colmat,sitemat)
row=[1;1;1;1;1;1;1;1;1;1];
col=[1;3;3;3;5;5;7;9;9;11];
site=[1;1;2;3;1;3;3;1;4;1];
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
%experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131217\';
experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';

rawdir = [imagepath,experimentpath,'Raw\'];

name1='DAPI_'; %nuc
name2='FITC_';
name3='Cy3_';
name4='Cy5_';

tempshot=[num2str(row(1)),'_',num2str(col(1)),'_',num2str(site(1))];
tempimg=single(imread([rawdir,tempshot,'_',name1,'stain.tif']));
[height,width]=size(tempimg);
raw1avg=zeros(height,width);
raw2avg=zeros(height,width);
raw3avg=zeros(height,width);
raw4avg=zeros(height,width);

numimg=numel(row);
filtrad=6;
for i=1:numimg
    shot=[num2str(row(i)),'_',num2str(col(i)),'_',num2str(site(i))];
    raw1=imfilter(single(imread([rawdir,shot,'_',name1,'stain.tif'])),fspecial('disk',filtrad),'symmetric');
    raw1avg=raw1avg+raw1;
    raw2=imfilter(single(imread([rawdir,shot,'_',name2,'stain.tif'])),fspecial('disk',filtrad),'symmetric');
    raw2avg=raw2avg+raw2;
    raw3=imfilter(single(imread([rawdir,shot,'_',name3,'stain.tif'])),fspecial('disk',filtrad),'symmetric');
    raw3avg=raw3avg+raw3;
    raw4=imfilter(single(imread([rawdir,shot,'_',name4,'stain.tif'])),fspecial('disk',filtrad),'symmetric');
    raw4avg=raw4avg+raw4;
end
raw1avg=raw1avg/numimg;
raw2avg=raw2avg/numimg;
raw3avg=raw3avg/numimg;
raw4avg=raw4avg/numimg;

bg1file=[rawdir,name1,'bgtemplate.tif'];
bg2file=[rawdir,name2,'bgtemplate.tif'];
bg3file=[rawdir,name3,'bgtemplate.tif'];
bg4file=[rawdir,name4,'bgtemplate.tif'];

fun=@(block_struct) nanmedian(block_struct.data(:));  %function handle for blockproc call
blocknum=10;
blockheight=ceil(height/blocknum); blockwidth=ceil(width/blocknum);
raw1block=blockproc(raw1avg,[blockheight blockwidth],fun);
raw1smooth=imresize(raw1block,[height width],'bicubic');
raw2block=blockproc(raw2avg,[blockheight blockwidth],fun);
raw2smooth=imresize(raw2block,[height width],'bicubic');
raw3block=blockproc(raw3avg,[blockheight blockwidth],fun);
raw3smooth=imresize(raw3block,[height width],'bicubic');
raw4block=blockproc(raw4avg,[blockheight blockwidth],fun);
raw4smooth=imresize(raw4block,[height width],'bicubic');

imwrite(uint16(raw1smooth),bg1file);
imwrite(uint16(raw2smooth),bg2file);
imwrite(uint16(raw3smooth),bg3file);
imwrite(uint16(raw4smooth),bg4file);
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