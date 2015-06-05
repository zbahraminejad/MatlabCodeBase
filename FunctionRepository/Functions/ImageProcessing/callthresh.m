function highthresh=callthresh(rawimage,iqrmult)
blurimg=imfilter(rawimage,fspecial('disk',10),'symmetric');
tiqr=prctile(blurimg(:),75)-prctile(blurimg(:),25);
tmean=mean(blurimg(:));
highthresh=tmean+iqrmult*tiqr;
