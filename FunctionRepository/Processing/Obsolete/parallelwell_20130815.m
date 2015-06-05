%close all;clear all; clear mex;
format short;
cd('H:\Documents\Timelapse\Code\Development\Processing')       % change to directory where program resides

%if ~matlabpool('size')
%    matlabpool('open')
%end

time1=tic;

%rows  = [3];
%cols  = [1 2 4 9 11];  %4 8 9 10
rows = [3 4];
cols = [3 5 6 7 8 10 12];
sites = [0 1];

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
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
    %%%%%% call any program with row,col,site %%%%%%%
    A_crop(row,col,site);
    %B_extractfeatures(row,col,site);
    %C_trackcells_replaceblurred(row,col,site);
    %D_tracesignals(row,col,site);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
toc(time1)
cd ..