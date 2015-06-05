function [x,y] = detectGlobalShift(im1,im2)
showImages=0; %flag for exploring images

%smooth images
smoothSize = 70;
smoothNorm = smoothSize^2;
kernel = ones(smoothSize,smoothSize)/smoothNorm;
smoothedIm1 = (conv2(double(im1),kernel,'valid'));
smoothedIm2 = (conv2(double(im2),kernel,'valid'));

%find cross correlation
z = abs(normxcorr2(smoothedIm1,smoothedIm2));

%find max correlation
[~, imax]=max(z(:));
[y,x] = ind2sub(size(z),imax);

%fix by size of image
x=(x-size(smoothedIm1,2));
y=(y-size(smoothedIm1,1));

if (showImages==1)
    figure;
    imagesc(z)
    title('cross corr results');
    T = maketform('affine', [1 0 0; 0 1 0; x y 1]);
    img2S = imtransform(smoothedIm1, T, ...
        'XData',[1 size(smoothedIm1,2)], 'YData',[1 size(smoothedIm1,1)]);

    img2 = imtransform(im1, T, ...
        'XData',[1 size(im1,2)], 'YData',[1 size(im1,1)]);

    figure ;
    subplot(3,2,1);
    imagesc((im1));
  
    title('im1');
    subplot(3,2,2);
    imagesc((im2));
  
    title('im2');
    subplot(3,2,3);
    imagesc(smoothedIm1);
    title('im1 - smoothed');
   
    subplot(3,2,4);
    imagesc(smoothedIm2);
    map = colormap();
    title('im2 - smoothed');
   
    
    
    subplot(3,2,5);
    imagesc((img2S));
    title('im1 - smoothed moved');
    
    subplot(3,2,6);
    imagesc((smoothedIm2));
    title('im2 - smoothed moved');
    
    figure;
    subplot(2,1,1);
    imshowpair(im1,im2);
    title('original overlay');
    subplot(2,1,2);
    imshowpair(img2,im2);
    title('Overlay after fix');
    
end




