function newmask=secondthresheach(image,blurradius,mask,boulderarea)
blur=imfilter(log(image),fspecial('disk',blurradius),'symmetric');
boulder_mask=bwareaopen(mask,boulderarea);
%bouldervals=image(boulder_mask);
[~,numboulders]=bwlabel(boulder_mask);
if numboulders==0
    newmask=mask;
    return;
end
boulder_info=regionprops(boulder_mask,'PixelIdxList');
%secondthresh_info=cell(numel(boulder_info),1);
secondthreshpixels=[];
for i=1:numboulders
    vals=blur(boulder_info(i).PixelIdxList);
    logvals=mat2gray(vals);
    secondthresh=graythresh(logvals);
    secondthresh=secondthresh*range(vals)+min(vals);
    secondthreshpixeach=boulder_info(i).PixelIdxList(vals>secondthresh);
    %secondthresh_info{i}=secondthreshpixels;
    secondthreshpixels=[secondthreshpixels;secondthreshpixeach];
end

[height,width]=size(image);
higher_mask=zeros(height,width);
higher_mask(secondthreshpixels)=1;
higher_mask(~boulder_mask)=0;
mask(boulder_mask)=0;
newmask=mask | higher_mask;
newmask=imfill(newmask,'holes');
end