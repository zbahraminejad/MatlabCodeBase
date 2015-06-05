function Average_ShadingImages(row,col,site)
%row=1;col=1;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadingpath= 'H:\Images\ShadingImages\20140402_DAPI_CFP_YFP_TxRed_Cy5\';

shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='DAPI_50ms_pink'; %nuc
name2='CFP_25ms_pink';
name3='YFP_25ms_yellow';
name4='TxRed_50ms_red';
name5='Cy5_100ms_red';
names={
    'DAPI_50ms_pink';
	'CFP_25ms_pink';
	'YFP_25ms_yellow';
	'TxRed_50ms_red';
	'Cy5_100ms_red';
    };
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
block1t=single(imread([shadingpath,'Raw\',shot,'_',name1,'.tif']));
sc1t=shadingcorrection_1([],block1,10);
block2t=single(imread([shadingpath,'Raw\',shot,'_',name2,'.tif']));
sc2t=shadingcorrection_1([],block2,10);
block3t=single(imread([shadingpath,'Raw\',shot,'_',name3,'.tif']));
sc3t=shadingcorrection_1([],block3,10);
block4t=single(imread([shadingpath,'Raw\',shot,'_',name4,'.tif']));
sc4t=shadingcorrection_1([],block4,10);
block5t=single(imread([shadingpath,'Raw\',shot,'_',name5,'.tif']));
sc5t=shadingcorrection_1([],block5,10);

%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'_stain.tif'])); raw1=raw1./sc1;
raw2=single(imread([rawdir,name2,'_stain.tif'])); raw2=raw2./sc2;
raw3=single(imread([rawdir,name3,'_stain.tif'])); raw3=raw3./sc3;
%raw4=single(imread([rawdir,name4,'_stain.tif'])); raw4=raw4./sc4;

if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
[nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %comment out for sparse YT analysis
%boulderarea=2000; nuc_mask=blobdetector_excludelargeandwarped(raw1,nucr,blobthreshold,debrisarea,boulderarea);
%%% remove objects that are highly eccentric %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% calculate background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dilated_mask=imdilate(foreground,strel('disk',nucr));
%nanmaskcyto=imdilate(foreground,strel('disk',2*nucr));
blur1=imfilter(raw1,fspecial('disk',nucr),'symmetric');
blur2=imfilter(raw2,fspecial('disk',nucr),'symmetric');
blur3=imfilter(raw3,fspecial('disk',nucr),'symmetric');
%blur4=imfilter(raw4,fspecial('disk',nucr),'symmetric');
iqrmult=4;
highfliers1=maskhighfliers_nonnuc(blur1,dilated_mask,nucr,iqrmult);
highfliers2=maskhighfliers_nonnuc(blur2,dilated_mask,nucr,iqrmult);
highfliers3=maskhighfliers_nonnuc(blur3,dilated_mask,nucr,iqrmult);
%highfliers4=maskhighfliers_nonnuc(blur4,dilated_mask,nucr,iqrmult);
fg1=dilated_mask | highfliers1;
fg2=dilated_mask | highfliers2;
fg3=dilated_mask | highfliers3;
%fg4=dilated_mask | highfliers4;
%fg=fg1;
fg=fg1 | fg2 | fg3;
%fg=fg1 | fg2 | fg3 | fg4;
blur1(fg)=NaN;
blur2(fg)=NaN;
blur3(fg)=NaN;
%blur4(fg)=NaN;
bg1=blocksmooth(blur1,10); %2nd Arg: blocknum
%bg1=blocksmooth_mode(blur1,10);
bg2=blocksmooth(blur2,10);
bg3=blocksmooth(blur3,10);
%bg4=blocksmooth(blur4,10);
%%% remove nucleoli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% feature_mask=smallblobdetector(log(raw2),12,-0.03,20,150);
% nuc_mask=logical(nuc_mask.*~feature_mask);
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real1=raw1-bg1;
real2=raw2-bg2;
real3=raw3-bg3;
%real4=raw4-bg4;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label=bwlabel(nuc_mask);
%smearmask=highfliers1;
smearmask=highfliers1 | highfliers2 | highfliers3;
%smearmask=highfliers1 | highfliers2 | highfliers3 | highfliers4;
smear_label=nuc_label.*smearmask;
nuc_smeared=unique(smear_label(:));
nuc_smeared(1)=[]; %remove 0
nuc_label(ismember(nuc_label,nuc_smeared))=0;
nuc_mask=nuc_label>0;
nuc_info=struct2cell(regionprops(nuc_mask,'Area','Centroid')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
if maskwrite
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=numel(nuc_area);
nuc_info=regionprops(nuc_mask,'PixelIdxList');
%ring_label=getcytoring(nuc_label);
%ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN;
sig1=nanvec;
sig2=nanvec; sig2_ring=nanvec;
sig3=nanvec; sig3_ring=nanvec;
sig4=nanvec;
for cc=1:numcells
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=mean(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=mean(real3(nuc_info(cc).PixelIdxList));
    %sig4(cc)=mean(real4(nuc_info(cc).PixelIdxList));
%     ringall=real2(ring_info(cc).PixelIdxList);
%      sig2_ring(cc)=mean(ringall(ringall>=median(ringall)));
%     ringall=real3(ring_info(cc).PixelIdxList);
%      sig3_ring(cc)=mean(ringall(ringall>=median(ringall)));
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
% IF_bleedthrough=bleedthroughrate(1)+RFP_sig*bleedthroughrate(2);
% IF_sig_med=IF_sig_med-IF_bleedthrough;
% IF_sig_max=IF_sig_max-IF_bleedthrough;
% IF_sig_mean=IF_sig_mean-IF_bleedthrough;
%%% correct DAPI for bleedthrough (use prior frame) %%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_CFPtoDAPI.mat'],'bleedthroughrate');
% nuc_sig_prev=tracedata(:,totalframes,4)./tracedata(:,totalframes,3); %mean CFP intensity
% DAPI_bleedthrough=bleedthroughrate(1)+nuc_sig_prev*bleedthroughrate(2);
% DAPI_sig=DAPI_sig-DAPI_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4,sig2_ring,sig3_ring];
save([datadir,shot,'.mat'],'IFdata');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}