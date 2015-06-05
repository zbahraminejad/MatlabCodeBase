%function C_trackcells(row,col,site)
cd ..\Functions; %change directory for function calls
%path='h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_10x_Steve\';
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = [path,'Data\'];
cpdir = [path,'CroppedProcessed\'];
nucr=8;

shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,'wellsss_', shot, '_MC.mat'],'wellsss');
rawimage = single(imread([cpdir,shot,'_nuc_1.tif']));    %open sample to get image size
dims = size(rawimage);

wellsp=relocateD(wellsss);  %links one object to the same object in the next frame
[newsp,trk_rc,corrupted,leaveout]=modifysp(wellsp,nucr,dims); %correction for merge and splitting and some gap filling
[bestsp,best_rc,corruptlist,leaveoutlist]=rmduplicates(newsp,trk_rc,corrupted,nucr,leaveout);
save([datadir,shot,'_alldata_MC.mat'],'wellsss','wellsp','newsp','trk_rc','bestsp','best_rc','corruptlist','leaveoutlist')
cd ..\Processing; %return to this directory