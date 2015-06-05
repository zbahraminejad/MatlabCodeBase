cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
cpdir = [path,'CroppedProcessed\'];
timepoint = '24hr';

rows = [1:8];
cols = [10];  %4:MK2206 10:DMSO
sites = [1];

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
cellcountmatrix=zeros(numrows,numcols,numsites);

nucroot = '_nuc_';
nucr=12;
minnucarea=round(pi*(nucr/4)^2);
frame=160;   %24hr=160; 0hr=40;

for rowidx=1:numrows
    for colidx=1:numcols
        for siteidx=1:numsites
time1=tic;
row=rows(rowidx); col=cols(colidx); site=sites(siteidx);
fprintf('row %0.0f, col %0.0f, site %0.0f\n',row,col,site);
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
DAs_or=single(imread([cpdir,shot,nucroot,num2str(frame),'.tif']));

%%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_bs=bgsub(log(DAs_or),10*nucr,0.05);
DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector
DAs_pad=bwareaopen(DAs_pad,minnucarea);
[~,numcells]=bwlabel(DAs_pad);
cellcountmatrix(rowidx,colidx,siteidx)=numcells;
toc(time1)
        end
    end
end
save([path,'cellcountmatrix_DMSO',timepoint,'_temp.mat'],'cellcountmatrix');
cd ..\Processing; %return to this directory