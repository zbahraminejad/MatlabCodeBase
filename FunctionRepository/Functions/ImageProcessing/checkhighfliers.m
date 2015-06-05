function [highflierflag,highthresh,blurimg]=callthresh(rawimage,iqrmult)
blurimg=imfilter(rawimage,fspecial('disk',50),'symmetric');
tiqr=prctile(blurimg(:),75)-prctile(blurimg(:),25);
tmean=mean(blurimg(:));
highthresh=tmean+iqrmult*tiqr;
highfliermask=blurimg>highthresh;
highflierflag=sum(highfliermask(:))>0;
