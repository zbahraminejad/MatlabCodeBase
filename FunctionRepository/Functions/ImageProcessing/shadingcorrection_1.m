function shadecorrected=shadingcorrection_1(emptyimg,blockimg,filtrad)
if isempty(emptyimg)
    emptyimg=zeros(size(blockimg));
end
shadecorrected=blockimg-emptyimg;
shadecorrected=imfilter(shadecorrected,fspecial('disk',filtrad),'symmetric');
maxillumination=prctile(shadecorrected(:),99);
shadecorrected=shadecorrected/maxillumination;