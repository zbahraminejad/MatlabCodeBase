%function C_trackcells_replaceblurred(row,col,site)
timetotal=tic;
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = [path,'Data\'];
cpdir = [path,'CroppedProcessed\'];
nucr=12;
nucsizslot=4;

shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,'wellsss_', shot, '.mat'],'wellsss');
rawimage = single(imread([cpdir,shot,'_nuc_1.tif']));    %open sample to get image size
dims = size(rawimage);

fmax=size(wellsss,3);
[wellsss,blurframes]=removeblurs(wellsss,nucsizslot,fmax);   %detect and remove blurs
wellsp=relocateD(wellsss);  %links one object to the same object in the next frame
[newsp,trk_rc,corrupted,leaveout]=modifysp_double(wellsp,nucr,dims); %correction for merge and splitting and some gap filling
[bestsp,best_rc,corruptlist,leaveoutlist]=rmduplicates(newsp,trk_rc,corrupted,nucr,leaveout);
[bestsp,best_rc]=replaceblurs(bestsp,best_rc,blurframes,fmax);

save([datadir,shot,'_alldata.mat'],'wellsss','wellsp','newsp','trk_rc','bestsp','best_rc','corruptlist','leaveoutlist','blurframes')
cd ..\Processing; %return to this directory
toc(timetotal)