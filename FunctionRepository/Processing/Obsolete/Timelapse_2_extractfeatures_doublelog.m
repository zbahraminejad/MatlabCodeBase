%function Timelapse_2_extractfeatures(row,col,site)
row='D';col='05';site='3';
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
datadir = ([path,'Data\']);
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
rawdir = [path,'Raw\',shot,'\'];
cpdir = [path,'Processed\',shot];
if ~exist(cpdir,'dir')
    mkdir(cpdir);
end

SF=229;EF=230; %Sabrina20x:208 Steve20x:218 Steve10x:110 Steve&Sabrina:240
initF=SF;   %the first intended frame (only matters if I have to restart
%OldAxon: MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%IX-Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12;
DAname = 'CFP_'; %nuc
DAedgename = 'nucedge&ring_';
REname = 'YFP_'; %DHB
CEname = 'TexasRed_'; %p21

debrisarea=pi*(nucr/4)^2;

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
continuation=0; restartframe=76;
if continuation==1
    load([datadir,'wellsss_',shot,'_restart'],'wellsss');
    wellsssrestart=cell(1,1,EF-SF+1);
    wellsssrestart(1:restartframe-1)=wellsss(1:restartframe-1);
    wellsss=wellsssrestart;
    SF=restartframe;
else
    wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)  
end
tempread=imread([rawdir,DAname,num2str(1),'.tif']);
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
    
    DAs_or=single(imread([rawdir,DAfile]));
    REs_or=single(imread([rawdir,REfile]));
    CEs_or=single(imread([rawdir,CEfile]));
    %%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    DAs_bs=bgsub(log(DAs_or),10*nucr,0.05);   %better for nuclear mask?
    REs_bs=bgsub_MC(log(REs_or),blockheight,blockwidth);
    CEs_bs=bgsub(log(CEs_or),10*nucr,0.05);  %local would be better--change
    %%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_pad=getnucmask(DAs_bs,nucr);
    %DAs_pad=getnucmask_onlybig(DAs_bs,nucr,2);
    
    %DAs_pad=getnucmask_histsweep_onlybig(DAs_bs,nucr,1,1);
    %[DAs_da,realnuc_la,finalcytoring]=simplercytoring(DAs_pad,nucr);
    %fcr_da=regionprops(finalcytoring,'PixelIdxList');
    
    %blob filter
    timecheck=tic;
    sigma=0.5*nucr/sqrt(2);
    h=sigma^2*fspecial('log',[nucr*5 nucr*5],sigma);
    nuc_log=imfilter(DAs_or,h,'symmetric');
    logmask=nuc_log<0;
    logmask=imopen(logmask,strel('disk',nucr/2,0));
    %imshow(logmask);
    
    sigma=0.1*nucr/sqrt(2);
    h=sigma^2*fspecial('log',[nucr*5 nucr*5],sigma);
    nuc_log=imfilter(DAs_or,h,'symmetric');
    logmaskfine=nuc_log<0;
    logmaskfine=~bwmorph(~logmaskfine,'diag');      %break connections
    logmaskfine=~bwmorph(~logmaskfine,'bridge');
    logmaskfine=bwareaopen(logmaskfine,round(pi*(nucr/3)^2));
    %imshow(logmaskfine);
    
    overlap=logmask_label.*logmaskfine;
    overlap=unique(overlap);
    overlap(1)=[]; %remove 0
    logmaskfinal=ismember(logmask_label,overlap);
    logmaskfinal=imfill(logmaskfinal,'holes');
    %imshow(logmaskfinal);
    
    %{
    tempframe=imadjust(mat2gray(DAs_bs));
    tempframe(:,:,2)=bwmorph(logmaskfinal,'remove');
    tempframe(:,:,3)=0;
    figure,imshow(tempframe);
    toc(timecheck)
    %}
    
    DAs_pad=detectdeflections(logmaskfinal,nucr,1);
    
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
    %extractmask=bwmorph(realnuc_la,'remove') + logical(finalcytoring);
    extractmask=bwmorph(realnuc_la,'remove');
    %extractmask=logical(finalcytoring);
    imwrite(uint16(extractmask),[cpdir,'\',NEfile]);
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
        
        %ring(cc,1)=median(REs_bs(fcr_da(cc).PixelIdxList));
        allringpixels=REs_bs(fcr_da(cc).PixelIdxList);
        topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
        ring(cc,1)=mean(topringpixels);
    end
    wellsss{f-initF+1}=[XX,YY,DD,AC,RR,CCC,ring]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame
    toc(timeframe)
end
save([datadir,'wellsss_', shot, '.mat'],'wellsss');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
%extractmask=bwmorph(DAs_pad,'remove');
tempframe=imadjust(mat2gray(DAs_bs));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
cd([codepath,'Processing\']); %return to this directory