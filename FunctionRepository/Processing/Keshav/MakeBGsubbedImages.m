function MakeBGsubbedImages
row=11;col=2;site=4;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadingpath='D:\Images\ShadingImages\20140410 DCYTC 20x\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];

rawdir = ['D:\Documents\Projects\Keshav\',shot,'_'];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='Hoechst'; %nuc
name2='DHB';
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG_bin2','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([shadingpath,'DAPI_bin2','.mat'],'shadingcorrection'); sc1=shadingcorrection;
load([shadingpath,'YFP_bin2','.mat'],'shadingcorrection'); sc2=shadingcorrection;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif'])); raw1=(raw1-bgcmos)./sc1;
raw2=single(imread([rawdir,name2,'_stain.tif'])); raw2=(raw2-bgcmos)./sc2;
%%% calculate background: min %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real1=bgsubmin(raw1,10);
real2=bgsubmin(raw2,10);
%%% write images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(uint16(real1),[rawdir,name1,'_bgsubbed.tif']);
imwrite(uint16(real2),[rawdir,name2,'_bgsubbed.tif']);