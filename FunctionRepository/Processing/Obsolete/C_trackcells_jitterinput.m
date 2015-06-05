%function C_trackcells_directory(row,col,site)
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

%row=3;col=3;site=1;
row='C';col='05';site='1';
%shot=[row,'_',col,'_',site];
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

%path='h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_10x_Steve\';
%path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%path = 'h:\Documents\Timelapse\Timescape\Steve_siCyclin\';
%path = 'h:\Documents\Timelapse\Timescape\20130207_Sabrina-MCF10A-p21KO\';
%path = 'h:\Documents\Timelapse\Timescape\Heewon_EKAR&DHB\';
%path = 'h:\Documents\Timelapse\Timescape\20130424_Panel_DHB-Gem_40hr\';
%path = 'H:\Documents\Timelapse\Timescape\20130514_MCF10A-p21KO_24hr_6min\';
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';

datadir = [path,'Data_Temp\'];
directory = 1;
if directory==1
    cpdir = [path,'Processed\',shot,'\'];
else
    cpdir = [path,'Processed\'];
end

nucr=12;
nucsizslot=4;

timetotal=tic;
load([datadir,'wellsss_', shot, '.mat'],'wellsss');
load([datadir,'jitter_', shot, '.mat'],'x','y');
rawimage = single(imread([cpdir,'nucedge&ring_1.tif']));    %open sample to get image size
%rawimage = single(imread([cpdir,shot,'_nuc_1.tif']));    %open sample to get image size

dims = size(rawimage);
fmax=size(wellsss,3);
wells_jittered=incorporatejitters(wellsss,x,y);
%[wells_noblurframes,blurframes]=removeblurs(wells_jittered,nucsizslot,fmax);   %detect and remove blurs
%if there were blurry frames, mask jitter calc has to take place afterwards
wellsp=relocateD(wells_noblurframes);  %links one object to the same object in the next frame
[newsp,trk_rc,corrupted,leaveout]=modifysp_double(wellsp,nucr,dims); %correction for merge and splitting and some gap filling
[bestsp_nbf,best_rc_nbf,corruptlist,leaveoutlist]=rmduplicates(newsp,trk_rc,corrupted,nucr,leaveout);
[bestsp,best_rc]=replaceblurs(bestsp_nbf,best_rc_nbf,blurframes,fmax);

save([datadir,shot,'_alldata.mat'],'bestsp','best_rc','corruptlist','leaveoutlist','blurframes')
toc(timetotal)
cd([codepath,'Processing\']); %return to this directory