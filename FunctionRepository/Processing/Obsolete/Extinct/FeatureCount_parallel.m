function [XX,YY,DD,AC,RR,CCC,avgringyfp] = FeatureCount_parallel(DAs_da, DAs_bs, REs_bs, CEs_bs, ringxypos)
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
end