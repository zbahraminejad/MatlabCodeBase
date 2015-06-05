cd('H:\Documents\Timelapse\Code\Development\Processing')       % change to directory where program resides
path = 'h:\Documents\Timelapse\Timescape\2013-04-18_DoseResponsePanel\';

time1=tic;

rows = [1:8];
cols = [1:12];
sites = [0:3];
timepoint = '0hr';
timepath = [path,timepoint,'\'];

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
cellcountmatrix=zeros(numrows,numcols,numsites);
shots=numrows*numcols*numsites;

parfor shot=1:shots
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
    fprintf('%0.0f_%0.0f_%0.0f\n',row,col,site);
    cellcountmatrix(rowidx,colidx,siteidx)=cellcount(row,col,site,timepath);    %get cell count
end
save([path,'cellcountmatrix_',timepoint,'.mat'],'cellcountmatrix');
toc(time1)
cd ..