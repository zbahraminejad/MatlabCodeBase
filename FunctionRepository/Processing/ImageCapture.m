function ImageCapture(row,col,site)
row=6;col=12;site=2;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
%experimentpath='2014-01-23_siCycA_CDK4i\';
%experimentpath='2014-01-13_Ab_Titrations\20140121_p21CycA\';
experimentpath='2014-01-13_Ab_Titrations\20140126_CycECycD1_p27p53_p57CycD2\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
reportdir = ([projectpath,experimentpath,'Report\']);
separatedirectories=0;
if separatedirectories
    rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
else
    rawdir = [imagepath,experimentpath,'Raw\',shot,'_'];
end
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='Hoechst'; %nuc
%name2='CycE';
name2='CycD1';
savename='CycD1_1to256000';
savefile=[reportdir,savename,'.tif'];
moviebin=1;
if moviebin==1
    nucr=12; %MCF-10A:12 YT+:8
    debrisarea=200; %MCF-10A:200 YT+:300
    boulderarea=10*debrisarea;
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
blobthreshold=-0.03; %YT+:-0.03
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif']));
raw2=single(imread([rawdir,name2,'_stain.tif']));

%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
[nuc_mask,~]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
%nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%boulderarea=2000; nuc_mask=blobdetector_excludelargeandwarped(raw1,nucr,blobthreshold,debrisarea,boulderarea);
%%% remove objects that are highly eccentric %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% view whole image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw2));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
%imshow(tempframe);
%%% save region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
miny=950; maxy=miny+350;
minx=1300; maxx=minx+350;
tempframe=tempframe(miny:maxy,minx:maxx,:);
imwrite(tempframe,savefile);