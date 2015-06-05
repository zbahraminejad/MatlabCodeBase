%function Timelapse_2_trackcells(row,col,site)
row='D';col='05';site='4';
%row='C';col='05';site='1';

projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath = '2013-06-07_p21_cy2_deletions\Experiment_20130715\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
datadir = ([projectpath,experimentpath,'Data_Test\']);
maskdir = [imagepath,experimentpath,'Mask_Test\',shot,'\'];

nucr=12;
immunostain=0;
stringIF='';IFstring='';
if immunostain==1
    stringIF='_IF';
    IFstring='IF_';
end

timetotal=tic;
load([datadir,'wellsss_',IFstring,shot,'.mat'],'wellsss','badframes','x','y'); %if IF, already added
fmax=size(wellsss,3);
xnobadframes=x; ynobadframes=y;
badframes=find(badframes==1);
wellsss(badframes)=[]; %remove bad frames from wellsss
rawimage = single(imread([maskdir,'nucedge_1.tif']));    %open sample to get image size

dims = size(rawimage);
wells_jittered=incorporatejitters(wellsss,xnobadframes,ynobadframes,dims);
wellsp=relocateD(wells_jittered);  %links one object to the same object in the next frame
[newsp,trk_rc,corrupted,leaveout]=modifysp_double(wellsp,nucr,dims); %correction for merge and splitting and some gap filling
[bestsp_ni,best_rc_ni,corruptlist,leaveoutlist]=rmduplicates(newsp,trk_rc,corrupted,nucr,leaveout);
[bestsp,best_rc,x,y]=interpolateframes(bestsp_ni,best_rc_ni,xnobadframes,ynobadframes,badframes,fmax);
save([datadir,shot,'_alldata',stringIF,'_nucmean.mat'],'bestsp','best_rc','corruptlist','leaveoutlist','badframes','x','y');
toc(timetotal);