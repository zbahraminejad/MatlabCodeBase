%close all;clear all; clear mex;
format short;
cd('H:\Documents\Timelapse\Code\Development\Processing')       % change to directory where program resides

%if ~matlabpool('size')
%    matlabpool('open')
%end

time1=tic;

rows = [1 2];
cols = [4 5 6];

numrows=length(rows);
numcols=length(cols);
shots=numrows*numcols;
parfor shot=1:shots
    colidx=mod(shot,numcols);
    if colidx==0
        colidx=numcols;
    end
    col=cols(colidx);
    rowidx=ceil(shot/numcols);
    row=rows(rowidx);
    fprintf('%0.0f_%0.0f\n',row,col);
    %%%%%% call any program with row,col,site %%%%%%%
    %A_crop_nosite(row,col);
    %B_extractfeatures(row,col);
    C_trackcells_replaceblurred(row,col);
    %D_tracesignals(row,col,site);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
toc(time1)
cd ..