function Immunostain_2(rowmat,colmat,sitemat)
row=3;col=2;site=2;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131217\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
rawdir = [imagepath,experimentpath,'Raw\'];

name1='DAPI_'; %nuc
name2='FITC_';
name3='Cy3_';
name4='Cy5_';

raw1=single(imread([rawdir,shot,'_',name1,'stain.tif']));
raw2=single(imread([rawdir,shot,'_',name2,'stain.tif']));
raw3=single(imread([rawdir,shot,'_',name3,'stain.tif']));
raw4=single(imread([rawdir,shot,'_',name4,'stain.tif']));
bg1file=[rawdir,name1,'bgtemplate.tif'];
bg2file=[rawdir,name2,'bgtemplate.tif'];
bg3file=[rawdir,name3,'bgtemplate.tif'];
bg4file=[rawdir,name4,'bgtemplate.tif'];

fun=@(block_struct) nanmedian(block_struct.data(:));  %function handle for blockproc call
[height,width]=size(raw1);
blocknum=10;
blockheight=ceil(height/blocknum); blockwidth=ceil(width/blocknum);
raw1block=blockproc(raw1,[blockheight blockwidth],fun);
raw1smooth=imresize(raw1block,size(raw1),'bicubic');
raw2block=blockproc(raw2,[blockheight blockwidth],fun);
raw2smooth=imresize(raw2block,size(raw2),'bicubic');
raw3block=blockproc(raw3,[blockheight blockwidth],fun);
raw3smooth=imresize(raw3block,size(raw3),'bicubic');
raw4block=blockproc(raw4,[blockheight blockwidth],fun);
raw4smooth=imresize(raw4block,size(raw4),'bicubic');

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