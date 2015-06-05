function ImageCorr=bgsub_MC(ImageOld,Radius,fra)
% fra: fraction; if zero, will be generated automatically

%% flatten the image
ImageNew=ImageOld-getbg_simple(ImageOld,Radius,fra);

%% get background value
%[a,b,bg]=ThreshImage(ImageNew);
[a,b,bg]=ThreshImage_MC(ImageNew,0);

%% subtract background value
ImageCorr=ImageNew-bg;
%ImageCorr(ImageCorr<0)=0;