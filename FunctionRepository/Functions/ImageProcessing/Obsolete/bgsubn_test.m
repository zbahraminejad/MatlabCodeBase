function [ImageCorr,Imagebg]=bgsubn_test(ImageOld,Imagemask,Radius)

%% preparation
ImageOld=single(ImageOld);
Imagemask=single(Imagemask);

%% flatten the image
Imagebg=getbg(ImageOld,~Imagemask,Radius);
ImageNew=ImageOld-Imagebg;

%% get background value
%[~,~,bg]=ThreshImage_smbg(ImageNew);
bg=0;

%% subtract background value
ImageCorr=ImageNew-bg;
ImageCorr(ImageCorr<0)=0;