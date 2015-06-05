imagepath='G:\Michael\';
experimentpath='20150401-CC-Diff\';
datadir=([imagepath,experimentpath,'Data\']);

rows=2:7;
cols=3;
sites=1:2;
numrows = length(rows);
numcols = length(cols);
numsites = length(sites);
shots = numrows*numcols*numsites;
nucr = 8;
motheroption = 0;
daughteroption = 0;
IFoption = 0;
figure
for shot = 1:shots
    
siteidx=mod(shot,numsites);
if siteidx==0
    siteidx=numsites;
end
site=sites(siteidx);
colidx=mod(ceil(shot/numsites),numcols);
if colidx==0
    colidx=numcols;
end
col=cols(colidx);
rowidx=ceil(shot/(numcols*numsites));
row=rows(rowidx);

shotID = [num2str(row),'_', num2str(col), '_', num2str(site)];
wellname = nameandsite(shotID);

rawdir=[imagepath,experimentpath,'Raw\',wellname,shotID,'_'];
datafile =[datadir,'tracedata_',shotID,'_nolink','.mat'];
load(datafile,'tracedata','genealogy','jitters');

[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother]=gathertracedata_mz_1(datadir,datafile,shot,motheroption,daughteroption,IFoption);

tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
subtightplot(3,4,shot),
reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
ylim([0 100]);

end