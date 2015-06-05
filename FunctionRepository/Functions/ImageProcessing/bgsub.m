function ImageCorr=bgsub(ImageOld,Radius,fra)
% fra: fraction; if zero, will be generated automatically

%% flatten the image
ImageNew=ImageOld-getbg_simple(ImageOld,Radius,fra);

%% get background value
[a,b,bg]=ThreshImage(ImageNew);

%% subtract background value
ImageCorr=ImageNew-bg;
%ImageCorr(ImageCorr<0)=0;    %This used to be commented out by FCT.  Uncomment for MC.