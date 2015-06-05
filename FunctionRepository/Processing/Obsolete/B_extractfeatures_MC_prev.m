%function B_extractfeatures_MC(row,col,site)
timetotal=tic;
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_10x_Steve\';  %folder containing the movies folders [CHANGE]
%path = 'h:\Documents\Timelapse\Timescape\20120807_12drugs\';
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
SF=1;EF=240; %Sabrina20x:208 Steve20x:218 Steve10x:110 Steve&Sabrina:240
EF=100;
initF=SF;   %the first intended frame (only matters if I have to restart
%OldAxon: MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%IX-Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12;
DAname = '_nuc_';
DAedgename = '_nucedge&ring_';
REname = '_hDHB_';
CEname = '_cdt1_';  %Pilot:'cdt1' Steve20x:'geminin' Steve10x:'geminin' Steve&Sabrina:'cdt1'
minnucarea=pi*(nucr/3)^2;
continuation=0; restartframe=76;

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
if continuation==1
    load([datadir,'wellsss_',shot,'_restart'],'wellsss');
    wellsssrestart=cell(1,1,EF-SF+1);
    wellsssrestart(1:restartframe-1)=wellsss(1:restartframe-1);
    wellsss=wellsssrestart;
    SF=restartframe;
else
    wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)  
end
tempread=imread([cpdir,shot,DAname,num2str(1),'.tif']);
[height,width]=size(tempread);
blocknum=3;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);

for f=SF:EF
    fprintf('frame %0.0f\n',f);
    timeframe=tic;

    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAfile=[shot,DAname,num2str(f),'.tif'];
    NEfile=[shot,DAedgename,num2str(f),'.tif'];
    REfile=[shot,REname,num2str(f),'.tif'];
    CEfile=[shot,CEname,num2str(f),'.tif'];  
    
    DAs_or=single(imread([cpdir,DAfile]));
    REs_or=single(imread([cpdir,REfile]));
    CEs_or=single(imread([cpdir,CEfile]));
    %%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    REs_bs=bgsub_MC(REs_or,blockheight,blockwidth);
    CEs_bs=bgsub_MC(CEs_or,blockheight,blockwidth);
    %%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector
    [DAs_da,finalcytoring]=buildcytoring(DAs_pad,REs_bs,nucr);
    fcr_da=regionprops(finalcytoring,'PixelIdxList');
    %%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcells=size(DAs_da,1);
    nucexclude=zeros(numcells,1);
    for k=1:numcells
        if DAs_da(k).Area < minnucarea
            nucexclude(k)=1;
        end
    end
    nucexclude=find(nucexclude);
    DAs_da(nucexclude)=[];
    fcr_da(nucexclude)=[];
    numcells=size(DAs_da,1);
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %realnuc_la(ismember(realnuc_la,nucexclude))=0;
    finalcytoring(ismember(finalcytoring,nucexclude))=0;
    extractmask=bwmorph(DAs_pad,'remove') + finalcytoring;
    %extractmask=finalcytoring;
    imwrite(uint16(extractmask),[cpdir,NEfile]);
    %%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tv=zeros(numcells,1);
    XX=tv; YY=tv; AC=tv;
    DD=tv;
    RR=tv;
    CCC=tv;
    ring=tv;
    for cc=1:numcells
        XX(cc,1)=DAs_da(cc).Centroid(1);  %x value of centroid
        YY(cc,1)=DAs_da(cc).Centroid(2);  %y value of centroid
        AC(cc,1)=DAs_da(cc).Area;
        
        DD(cc,1)=median(DAs_bs(DAs_da(cc).PixelIdxList));
        CCC(cc,1)=median(CEs_bs(DAs_da(cc).PixelIdxList));
        
        RR(cc,1)=median(REs_bs(DAs_da(cc).PixelIdxList));
        ring(cc,1)=median(REs_bs(fcr_da(cc).PixelIdxList));
        
    end
    wellsss{f-initF+1}=[XX,YY,DD,AC,RR,CCC,ring]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame
    toc(timeframe)
end
save([datadir,'wellsss_', shot, '_HS_OP.mat'],'wellsss');
cd ..\Processing; %return to this directory
toc(timetotal);