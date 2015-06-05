function B_extractfeatures_directory_nocdt1(row,col,site)
%row='A';col='01';site='1';
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

shot=[row,'_',col,'_',site];

%path = 'h:\Documents\Timelapse\Timescape\20130424_Panel_DHB-Gem_40hr\';
path = 'H:\Documents\Timelapse\Timescape\20130514_MCF10A-p21KO_24hr_6min\';


datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\',shot,'\'];
SF=1;EF=205; %Sabrina20x:208 Steve20x:218 Steve10x:110 Steve&Sabrina:240
initF=SF;   %the first intended frame (only matters if I have to restart
%OldAxon: MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%IX-Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12;
nucname = 'nuc_';
nucedgename = 'nucedge&ring_';
DHBname = 'DHB_';
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
tempread=imread([cpdir,nucname,num2str(1),'.tif']);
[height,width]=size(tempread);
blocknum=3;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);

timetotal=tic;
for f=SF:EF
    fprintf('frame %0.0f\n',f);
    timeframe=tic;

    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucfile=[nucname,num2str(f),'.tif'];
    edgefile=[nucedgename,num2str(f),'.tif'];
    DHBfile=[DHBname,num2str(f),'.tif'];
    
    nuc_or=single(imread([cpdir,nucfile]));
    DHB_or=single(imread([cpdir,DHBfile]));
    %%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    nuc_bs=bgsub(log(nuc_or),10*nucr,0.05);   %better for nuclear mask?
    DHB_bs=bgsub_MC(log(DHB_or),blockheight,blockwidth);
    %%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucmask=getnucmask_histsweep(nuc_bs,nucr);  %MC histogram sweep & concave detector
    [DHB_info,nuc_labeled,finalcytoring]=buildcytoring(nucmask,DHB_bs,nucr);
    ring_info=regionprops(finalcytoring,'PixelIdxList');
    %%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcells=size(DHB_info,1);
    numrings=size(ring_info,1);
    if numcells>numrings
        nucexclude=[numrings+1:numcells];
        DHB_info(nucexclude)=[];
    end
    nucexclude=zeros(numrings,1);
    for k=1:numrings
        if DHB_info(k).Area < debrisarea
            nucexclude(k)=1;
        end
    end
    nucexclude=find(nucexclude);
    DHB_info(nucexclude)=[];
    ring_info(nucexclude)=[];
    numcells=size(DHB_info,1);

    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_labeled(ismember(nuc_labeled,nucexclude))=0;
    finalcytoring(ismember(finalcytoring,nucexclude))=0;
    extractmask=bwmorph(nuc_labeled,'remove') + logical(finalcytoring);
    %extractmask=bwmorph(realnuc_la,'remove');
    imwrite(uint16(extractmask),[cpdir,edgefile]);
    %%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tv=zeros(numcells,1);
    XX=tv; YY=tv; AC=tv;
    nucdata=tv;
    DHBnucdata=tv;
    DHBringdata=tv;
    for cc=1:numcells
        XX(cc,1)=DHB_info(cc).Centroid(1);  %x value of centroid
        YY(cc,1)=DHB_info(cc).Centroid(2);  %y value of centroid
        AC(cc,1)=DHB_info(cc).Area;
        
        nucdata(cc,1)=median(nuc_bs(DHB_info(cc).PixelIdxList));
        
        DHBnucdata(cc,1)=median(DHB_bs(DHB_info(cc).PixelIdxList));
        DHBringdata(cc,1)=median(DHB_bs(ring_info(cc).PixelIdxList));
        
    end
    wellsss{f-initF+1}=[XX,YY,nucdata,AC,DHBnucdata,DHBringdata]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame
    save([datadir,'wellsss_', shot, '.mat'],'wellsss');
    toc(timeframe)
end
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
tempframe=imadjust(mat2gray(DHB_bs));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}
cd([codepath,'Processing\']); %return to this directory