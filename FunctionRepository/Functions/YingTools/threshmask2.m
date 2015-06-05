function [mask, blur]=threshmask2(image, size)
blur=imfilter(image,fspecial('disk',size),'symmetric');
normlog=mat2gray(log(blur));
thresh=graythresh(normlog);
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end