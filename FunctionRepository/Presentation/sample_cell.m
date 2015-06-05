function sample_cell()
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
%%% set average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sample specification %%%%%%%%%%%
row=1;
col=11;
site=1;
track=214;
%%% Setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,shot,'_alldata_HS_L05'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
windowsize = 20;
windowhalf = windowsize/2;

%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[height,width]=size(tempread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=VideoWriter([datadir,shot,'_',track],'Uncompressed AVI');
M.FrameRate = 4;
open(M);

SF=best_rc(track,1);
EF=best_rc(track,3);
frames=SF:EF;

signal1 = -10000*ones(1,EF);
signal2 = -10000*ones(1,EF);

for f=SF:EF
    signal1(f) = bestsp{f}(track,7)./bestsp{f}(track,5);
    signal2(f) = bestsp{f}(track,6);
    
    %%% Read in pre-cropped images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    NEs_or=single(imread([cpdir,shot,'_nucedge&ring_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    CEs_or=single(imread([cpdir,shot,'_geminin_',num2str(f),'.tif']));
    
    %%%%%  Focus on track   %%%%%%%%%%%%%%
    cx=int16(bestsp{f}(track,2));
    cy=int16(bestsp{f}(track,1));
    miny=cy-(windowhalf-1)*nucr; maxy=cy+windowhalf*nucr;
    minx=cx-(windowhalf-1)*nucr; maxx=cx+windowhalf*nucr;
    if minx<1
        minx=1; maxx=1+windowsize*nucr;
    end
    if miny<1
        miny=1; maxy=1+windowsize*nucr;
    end
    if maxx>height
        maxx=height; minx=height-windowsize*nucr;
    end
    if maxy>width
        maxy=width; miny=width-windowsize*nucr;
    end
    DAs_or=DAs_or(minx:maxx,miny:maxy);
    NEs_or=NEs_or(minx:maxx,miny:maxy);
    REs_or=REs_or(minx:maxx,miny:maxy);
    CEs_or=CEs_or(minx:maxx,miny:maxy);

    %%% Add frame to movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempframe=imadjust(mat2gray(REs_or));
    tempframe(:,:,2)=imadjust(mat2gray(NEs_or));
    tempframe(:,:,3)=0;
    %writeVideo(M,im2frame(tempframe));
end

%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(frames,signal1);

cd ..\Analysis; %return to this directory

end