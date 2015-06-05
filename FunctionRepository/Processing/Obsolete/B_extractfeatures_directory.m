function B_extractfeatures_directory(row,col,site)
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

%row='B';col='04';site='2';
shot=[row,'_',col,'_',site];
%shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_10x_Steve\';  %folder containing the movies folders [CHANGE]
%path = 'h:\Documents\Timelapse\Timescape\20120807_12drugs\';
%path = 'h:\Documents\Timelapse\Timescape\Steve_siCyclin\';
%path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%path = 'h:\Documents\Timelapse\Timescape\20130424_Panel_DHB-Gem_40hr\';
%path = 'h:\Documents\CDK4\20130715\';
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';

datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\',shot,'\'];
SF=1;EF=230; %Sabrina20x:208 Steve20x:218 Steve10x:110 Steve&Sabrina:240
initF=SF;   %the first intended frame (only matters if I have to restart
%OldAxon: MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%IX-Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12;
DAname = 'nuc_';
DAedgename = 'nucedge&ring_';
REname = 'DHB_';
CEname = 'p21_';  %Pilot:'cdt1' Steve20x:'geminin' Steve10x:'geminin' Steve&Sabrina:'cdt1'
debrisarea=pi*(nucr/4)^2;
continuation=0; restartframe=76;

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if continuation==1
    load([datadir,'wellsss_',shot,'_restart'],'wellsss');
    wellsssrestart=cell(1,1,EF-SF+1);
    wellsssrestart(1:restartframe-1)=wellsss(1:restartframe-1);
    wellsss=wellsssrestart;
    SF=restartframe;
else
    wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)  
end
tempread=imread([cpdir,DAname,num2str(1),'.tif']);
[height,width]=size(tempread);
blocknum=3;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);

timetotal=tic;
for f=SF:EF
    fprintf('frame %0.0f\n',f);
    timeframe=tic;

    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAfile=[DAname,num2str(f),'.tif'];
    NEfile=[DAedgename,num2str(f),'.tif'];
    REfile=[REname,num2str(f),'.tif'];
    CEfile=[CEname,num2str(f),'.tif'];  
    
    DAs_or=single(imread([cpdir,DAfile]));
    REs_or=single(imread([cpdir,REfile]));
    CEs_or=single(imread([cpdir,CEfile]));
    %%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    DAs_bs=bgsub(log(DAs_or),10*nucr,0.05);   %better for nuclear mask?
    REs_bs=bgsub_MC(log(REs_or),blockheight,blockwidth);
    CEs_bs=bgsub(log(CEs_or),10*nucr,0.05);  %local would be better--change
    %%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector
    DAs_pad=getnucmask(DAs_bs,nucr);
    [DAs_da,realnuc_la,finalcytoring]=buildcytoring(DAs_pad,REs_bs,nucr);
    fcr_da=regionprops(finalcytoring,'PixelIdxList');
    %%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcells=size(DAs_da,1);
    numrings=size(fcr_da,1);
    if numcells>numrings
        nucexclude=[numrings+1:numcells];
        DAs_da(nucexclude)=[];
    end
    nucexclude=zeros(numrings,1);
    for k=1:numrings
        if DAs_da(k).Area < debrisarea
            nucexclude(k)=1;
        end
    end
    nucexclude=find(nucexclude);
    DAs_da(nucexclude)=[];
    fcr_da(nucexclude)=[];
    numcells=size(DAs_da,1);

    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    realnuc_la(ismember(realnuc_la,nucexclude))=0;
    finalcytoring(ismember(finalcytoring,nucexclude))=0;
    extractmask=bwmorph(realnuc_la,'remove') + logical(finalcytoring);
    %extractmask=bwmorph(realnuc_la,'remove');
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
    save([datadir,'wellsss_', shot, '.mat'],'wellsss');
    toc(timeframe)
end
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(REs_bs));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}
cd([codepath,'Processing\']); %return to this directory