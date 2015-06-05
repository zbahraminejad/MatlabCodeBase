function Timelapse(row,col,site)
row=3;col=4;site=1;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadingpath='H:\Images\ShadingImages\20140410_DCYTC_20x\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
blurrad=10;

load([shadingpath,'BG_',num2str(site),'.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([shadingpath,'CFP_',num2str(site),'.mat'],'shadingcorrection'); sc1=shadingcorrection;
load([shadingpath,'YFP_',num2str(site),'.mat'],'shadingcorrection'); sc2=shadingcorrection;

%raw1=single(imread([rawdir,name1,num2str(f),'.tif']));
raw1=single(imread([shadingpath,'Raw\',shot,'_CFP_serum.tif']));
raw1=imfilter(raw1,fspecial('disk',blurrad),'symmetric');
raw1=raw1-bgcmos;
raw1=raw1./sc1;
raw2=single(imread([shadingpath,'Raw\',shot,'_YFP_serum.tif']));
raw2=imfilter(raw2,fspecial('disk',blurrad),'symmetric');
raw2=raw2-bgcmos;
raw2=raw2./sc2;
keyboard;
imagesc(raw2);colorbar;
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=imadjust(mat2gray(raw1));
tempframe(:,:,3)=0;
tempframe(:,:,1)=extractmask;
figure,imshow(tempframe);
%}