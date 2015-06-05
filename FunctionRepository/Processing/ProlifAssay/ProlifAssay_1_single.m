cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
path = 'h:\Documents\Timelapse\Timescape\Heewon_EKAR&DHB\';
cpdir = [path,'Raw\'];

%rows = [1:8];
%cols = [10];  %4:MK2206 10:DMSO
%sites = [1];
row=3; col=6; site=1;

nucroot = '_CFP_';
nucr=12;
minnucarea=round(pi*(nucr/4)^2);
frame=[40 160];   %0hr=40; 24hr=160;
numcells = zeros(2,1);

for idx=1:2
    shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
    DAs_or=single(imread([cpdir,shot,nucroot,num2str(frame(idx)),'.tif']));
    DAs_bs=bgsub(log(DAs_or),10*nucr,0.05);
    DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector
    DAs_pad=bwareaopen(DAs_pad,minnucarea);
    [~,numcells(idx)]=bwlabel(DAs_pad);
end

ndiff=numcells(2)-numcells(1);
foldchange=ndiff/numcells(1)+1;    %normalized difference

fprintf('0hr cell count = %0.0f\n',numcells(1));
fprintf('24hr cell count = %0.0f\n',numcells(2));
fprintf('fold change = %0.2f\n',foldchange);
cd ..\Processing; %return to this directory