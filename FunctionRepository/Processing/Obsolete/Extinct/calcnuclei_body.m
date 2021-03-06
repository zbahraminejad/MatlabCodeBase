function main_body(row,col,site,SF,EF,tempdir,rawdir,nucr)

%% setup tempwell
movieName=[num2str(row),'_', num2str(col), '_', num2str(site)];  %[CHANGE]
disp(movieName);
%cropamount=50;
wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)
M=avifile([tempdir,'\Data\',movieName, '.avi'],'compression','none','fps',4);

%% correcting jitters
x=zeros(1,EF); y=zeros(1,EF);

%%
for f=SF:EF
    disp(f) %print to the screen the frame you are on
    time2=tic;

    %% reading images
    DAs_or=single(imread([rawdir,'\',movieName, '_mRFP1_', num2str(f), '.tif']));%DA = dapi, reads image. 2D matrix of intensities for that frame

    %% correct multiple jitters
    DI1=imadjust(mat2gray(log(DAs_or(1:925, :))));  %only use the first 925 rows of pixels bc sometimes images have black shadows/lines at the bottom

    if f>SF                   
        [xx,yy]=CalcJitter(DI0,DI1);
        x(f)=x(f-1)+xx;  
        y(f)=y(f-1)+yy;
        disp([x(f), y(f)]);
    end
    DI0=DI1;
    toc(time2);
end

cropamountx=ceil(max(abs(x))+1);  %use maximum jitter to find amount to crop  
cropamounty=ceil(max(abs(y))+1);  %use maximum jitter to find amount to crop         

for f=SF:EF
    time2=tic;

    %% reading images
    DAs_or=single(imread([rawdir,'\',movieName, '_mRFP1_', num2str(f), '.tif']));%DA = dapi, reads image. 2D matrix of intensities for that frame
    REs_or=single(imread([rawdir,'\',movieName, '_EYFP_',  num2str(f), '.tif']));%RE = red, reads image
    CEs_or=single(imread([rawdir,'\',movieName, '_ECFP_',  num2str(f), '.tif']));%CE = cerulean, reads image


    %% cropping.  

    DAs_or=CropJitter(DAs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));  %75 is my guess at the biggest jitter in pixels in whole movie
    REs_or=CropJitter(REs_or,   cropamountx, cropamountx, cropamounty, cropamounty,   x(f), y(f));  %75 is my guess at the biggest jitter in pixels in whole movie
    CEs_or=CropJitter(CEs_or,  cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));  %75 is my guess at the biggest jitter in pixels in whole movie


    %% image processing
    DAs_bl=log(imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric')); %take log to decrease the variance of the signal.  use fspecial to create a disk-shaped filter to make a blurred (bl) image
    DAs_bs=bgsub(DAs_bl,10*nucr,0.05);  %background subtration on the blurred image. search 10x the radius, sort, find the value of the bottom 5th percentile. can change these #s
    REs_bl=(imfilter(REs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    REs_bs=bgsub(REs_bl,10*nucr,0.05);
    CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    CEs_bs=bgsub(CEs_bl,10*nucr,0.05);


    %% Make cropped avi movie to feed into Plot_Nuclei_2c_withIF

    tempframe=imadjust(mat2gray(DAs_bs));  %red; to make a movie incluing nuclear marker
    tempframe(:,:,2)=imadjust(mat2gray(REs_bs)); %green.  
    tempframe(:,:,3)=imadjust(mat2gray(CEs_bs)); %blue.  

    M=addframe(M,im2frame(tempframe));

    %% get data
    DAs_ma=getdapimask(DAs_bs,nucr);  %make the mask. this segmentation function is the engine of Feng Chiao's script
    DAs_la=bwlabel(DAs_ma);  %labels the objects with numbers

    dapimaskdilated=imdilate(DAs_la, strel('disk', 1, 8));  %dilate dapimask by 2 pixels; 8 means octagon   
    ring=imdilate(dapimaskdilated, strel('disk', 3, 8)) - dapimaskdilated;  %dilate by 2 pixels, and subtract to make a ring              
    DAs_da=regionprops(DAs_la,'Centroid','PixelIdxList','Area'); %finds the centroid,etc of each labeled object  %to test, type "DAs_da(1).Area"
    ringxypos=regionprops(ring, 'PixelIdxList');


    %%
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
      %% filter data.  filters out CFP that has leaked into the cytoplasm and is now being recognized as an object

 %%  use this if you don't need to filter the data 
      %wellsss{:,:,f-SF+1}=[XX,YY,DD,AC,RR,CCC,avgringyfp]; %matrix.  rows are cells.  columns are x coor, y coord, nuclear intensity, area of nucleus, reporter intensity. for this frame
      wellsss{:,:,f}=[XX,YY,DD,AC,RR,CCC,avgringyfp];
      toc(time2)
end

%%
save([tempdir,'\Data\wellsss_', movieName, '.mat'],'wellsss')
M=close(M);

%size(DAs_or)  %display the size of the final cropped movie, to enter into Plot_Nuclei function
end
