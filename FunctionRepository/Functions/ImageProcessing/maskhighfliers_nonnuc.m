function mask=callthresh(rawimage,nuc_mask,nucr,iqrmult)
%rawimage=single(imread(['H:\Images\2013-11-26_CycDp21SerumRelease\Experiment_20131216\Raw\2_1_1_DAPI_stain.tif']));
%blurimg=imfilter(rawimage,fspecial('disk',10),'symmetric');
blurimg=rawimage;
nonnucintensities=blurimg(nuc_mask==0);
blurimg(nuc_mask)=0;
tiqr=prctile(nonnucintensities(:),75)-prctile(nonnucintensities(:),25);
tmean=mean(nonnucintensities(:));
highthresh=tmean+iqrmult*tiqr;
mask=blurimg>highthresh;
mask=imdilate(mask,strel('disk',round(1.5*nucr)));
