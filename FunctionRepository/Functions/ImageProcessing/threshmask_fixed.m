function mask=threshmask(image,blurradius)
blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
% normlog=mat2gray(log(blur));
norm = mat2gray(image);
thresh=graythresh(norm);
mask=im2bw(norm,thresh);
mask = imopen(mask,strel('disk',6));
mask=imfill(mask,'holes');
end