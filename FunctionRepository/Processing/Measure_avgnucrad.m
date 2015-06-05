%function Measure_avgnucrad
row='A';col='02';site='1';

imagepath='H:\Images\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=7;
frame=220;
nucname='CFP_'; %nuc
%%% segment and count cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_raw=single(imread([rawdir,nucname,num2str(frame),'.tif']));
nuc_mask=blobdetector(log(nuc_raw),nucr,-0.03);
nuc_mask=segmentdeflections(nuc_mask,nucr,0);
%%% get area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=cell2mat(struct2cell(regionprops(nuc_mask,'Area')));
figure, hist(nuc_area,100), title('nuclear area');
%%% calculate average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radii=sqrt(nuc_area/pi);
figure, hist(radii,100), title('nuclear radii');
avgrad=round(mean(radii));
medrad=round(median(radii));
fprintf('average radius = %0.0f\n',avgrad);
fprintf('median radius = %0.0f\n',medrad);
%%% debugging: view images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}