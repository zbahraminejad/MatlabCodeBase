function Timelapse_3_trackcells(row,col,site)
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

%row=3;col=3;site=1;
%row='D';col='05';site='4';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\';
datadir = [path,'Data\'];
cpdir = [path,'Processed\',shot,'\'];

nucr=12;

timetotal=tic;
immunostain=1;
IFstring=''; stringIF='';
if immunostain==1
    IFstring='IF_';
    stringIF='_IF';
end
load([datadir,'wellsss_',IFstring,shot,'.mat'],'wellsss','badframes'); %if IF, already added
fmax=size(wellsss,3);
load([datadir,'jitter_',IFstring,shot,'.mat'],'x','y');
xnobadframes=x; ynobadframes=y;
badframes=find(badframes==1);
wellsss(badframes)=[]; %remove bad frames from wellsss
rawimage = single(imread([cpdir,'nucedge&ring_1.tif']));    %open sample to get image size

dims = size(rawimage);
wells_jittered=incorporatejitters(wellsss,xnobadframes,ynobadframes,dims);
wellsp=relocateD(wells_jittered);  %links one object to the same object in the next frame
[newsp,trk_rc,corrupted,leaveout]=modifysp_double(wellsp,nucr,dims); %correction for merge and splitting and some gap filling
[bestsp_ni,best_rc_ni,corruptlist,leaveoutlist]=rmduplicates(newsp,trk_rc,corrupted,nucr,leaveout);
[bestsp,best_rc,x,y]=interpolateframes(bestsp_ni,best_rc_ni,xnobadframes,ynobadframes,badframes,fmax);
save([datadir,shot,'_alldata',stringIF,'.mat'],'bestsp','best_rc','corruptlist','leaveoutlist','badframes')
save([datadir,'jitter_',IFstring,shot,'.mat'],'x','y');
toc(timetotal)
cd([codepath,'Processing\']); %return to this directory