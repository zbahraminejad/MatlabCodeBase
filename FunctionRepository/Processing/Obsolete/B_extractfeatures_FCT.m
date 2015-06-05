%function B_extractfeatures(row,col,site)
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
SF=1;EF=240; %Sabrina20x:208 Steve20x:218 Steve&Sabrina10x:240
SF=1;EF=100;
nucr=8; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
DAname = '_nuc_';
DAedgename = '_nucedge&ring_FCT_';
REname = '_hDHB_';
CEname = '_cdt1_';

%%% setup tempwell
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)    

for f=SF:EF
    time2=tic;

    %%% reading images
    DAfile=[shot,DAname,num2str(f),'.tif'];
    NEfile=[shot,DAedgename,num2str(f),'.tif'];
    REfile=[shot,REname,num2str(f),'.tif'];
    CEfile=[shot,CEname,num2str(f),'.tif'];  
    
    DAs_or=single(imread([cpdir,DAfile]));
    REs_or=single(imread([cpdir,REfile]));
    CEs_or=single(imread([cpdir,CEfile])); 

    %%% image processing
    DAs_bl=log(imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric')); %take log to decrease the variance of the signal.  use fspecial to create a disk-shaped filter to make a blurred (bl) image
    DAs_bs=bgsub(DAs_bl,10*nucr,0.05);  %background subtraction on the blurred image. search 10x the radius, sort, find the value of the bottom 5th percentile. can change these #s
    REs_bl=(imfilter(REs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    REs_bs=bgsub(REs_bl,10*nucr,0.05);
    CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    CEs_bs=bgsub(CEs_bl,10*nucr,0.05);

    %%% get data
    DAs_ma=getdapimask(DAs_bs,nucr);  %make the mask. this segmentation function is the engine of Feng Chiao's script
    %nuclearedges = bwmorph(DAs_ma,'remove');
    DAs_la=bwlabel(DAs_ma);  %labels the objects with numbers
    dapimaskdilated=imdilate(DAs_la, strel('disk', 1, 8));  %dilate dapimask by 2 pixels; 8 means octagon   
    ring=imdilate(dapimaskdilated, strel('disk', 3, 8)) - dapimaskdilated;  %dilate by 2 pixels, and subtract to make a ring              
    DAs_da=regionprops(DAs_la,'Centroid','PixelIdxList','Area'); %finds the centroid,etc of each labeled object  %to test, type "DAs_da(1).Area"
    ringxypos=regionprops(ring, 'PixelIdxList');        
    %imwrite(uint16(nuclearedges),[cpdir,NEfile]);
    
    %%%
    mz = size(DAs_da,1);
    XX=zeros(mz,1);YY=zeros(mz,1);  %define vectors which will hold the info from the structured array, DAs_da 
    AC=zeros(mz,1);%PP=zeros(size(DAs_da,1),1);DI=zeros(size(DAs_da,1),1);
    DD=zeros(mz,1);RR=zeros(mz,1); CCC=zeros(mz,1);  %FR=zeros(size(DAs_da,1),1); %DD is a vector of median dapi intensities for each object in log scale
    avgringyfp=zeros(mz,1);
    for cc=1:mz   %run a loop to put the info from the structured array in to the vectors; this is for each cell
        XX(cc,1)=DAs_da(cc).Centroid(1);  %x value of centroid
        YY(cc,1)=DAs_da(cc).Centroid(2);  %y value of centroid
        AC(cc,1)=DAs_da(cc).Area;
        DD(cc,1)=median(DAs_bs(DAs_da(cc).PixelIdxList));
        RR(cc,1)=median(REs_bs(DAs_da(cc).PixelIdxList));  %nuc yfp intensity
        CCC(cc,1)=median(CEs_bs(DAs_da(cc).PixelIdxList));  %nuc cfp intensity

        allringpixels=REs_bs(ringxypos(cc).PixelIdxList);
        topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
        avgringyfp(cc,1)=mean(topringpixels);  %take mean of pixels in the top 50th pctile
    end

    %%%  use this if you don't need to filter the data 
    wellsss{f-SF+1}=[XX,YY,DD,AC,RR,CCC,avgringyfp]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame
    toc(time2)
end
save([datadir,'wellsss_', shot, '.mat'],'wellsss')
cd ..\Processing; %return to this directory