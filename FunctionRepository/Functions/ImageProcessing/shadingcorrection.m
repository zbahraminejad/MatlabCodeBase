function shadecorrected=shadingcorrection(emptyimg,blockimg,filtrad)
shadecorrected=blockimg-emptyimg;
shadecorrected=imfilter(shadecorrected,fspecial('disk',filtrad),'symmetric');
maxillumination=prctile(shadecorrected(:),99);
shadecorrected=shadecorrected/maxillumination;