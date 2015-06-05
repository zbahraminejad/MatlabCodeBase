function tracenuclei(row,col,site)
path='h:\Documents\Timescape\20120807_12drugs\';
datadir = [path,'Data\'];
nucr=8;  %change this

movieName=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,'wellsss_', movieName, '.mat'],'wellsss');
wellsp=relocateD(wellsss);  %links one object to the same object in the next frame
[newsp,trk_rc]=modifysp(wellsp); %correction for merge and splitting and some gap filling
[bestsp,best_rc]=rmduplicates(newsp,trk_rc,nucr);
save([datadir,movieName,'_alldata.mat'],'wellsss','wellsp','newsp','trk_rc','bestsp','best_rc')

end


