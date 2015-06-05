rows = 2:7;
cols = [3:5];
sites = 1:5;
days = [1:4];
numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
numdays = length(days);
shots=numrows*numcols*numsites;
for day = 1:numdays
parfor shot = 1:shots
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
        
% Fixed_Cell_analysis(row,col,site)
Fixed_Cell_analysis_multiple_day(row,col,site,days(day));
% fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'\n']);

end
end