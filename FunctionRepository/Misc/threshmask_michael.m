function mask=threshmask(image)
blur=imfilter(image,fspecial('disk',3),'symmetric');
% blur=blur-imfilter(blur,fspecial('disk',30),'symmetric');

normlog=mat2gray(log(blur));
% normlog=mat2gray(blur);
thresh=graythresh(normlog);
mask=im2bw(normlog,thresh);
mask=imfill(mask,'holes');
end