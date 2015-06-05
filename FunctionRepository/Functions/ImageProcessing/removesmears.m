function blurimage=removesmears(rawimage,highthresh)
blurimage=imfilter(rawimage,fspecial('disk',50),'symmetric'); %originally 10
mask=blurimage>highthresh;
if max(mask(:))>0
    mask=imdilate(mask,strel('disk',50));
    blurimage(mask)=NaN;
end