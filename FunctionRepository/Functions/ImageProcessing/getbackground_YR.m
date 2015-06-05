function [YFP_bg,RFP_bg]=getbackground_YR(nuc_mask,centers,YFP_raw,RFP_raw,nucr,ratio,YFPprc,RFPprc)
nuc_mask=imresize(nuc_mask,ratio,'nearest');
YFP_raw=imresize(YFP_raw,ratio,'nearest');
RFP_raw=imresize(RFP_raw,ratio,'nearest');

offset=ratio*20*nucr;
regionlength=2*offset+1;
[height,width]=size(nuc_mask);
YFP_raw(nuc_mask>0)=NaN;
RFP_raw(nuc_mask>0)=NaN;

centers=centers(~isnan(centers(:,1)),:);
centers=centers*ratio;
numcells=size(centers,1);
x=centers(:,1); y=centers(:,2);

onesmatrix=ones(numcells,1);
minrow=round(y-offset); minrow=max(minrow,onesmatrix);
maxrow=round(y+offset); maxrow=min(maxrow,onesmatrix*height);
mincol=round(x-offset); mincol=max(mincol,onesmatrix);
maxcol=round(x+offset); maxcol=min(maxcol,onesmatrix*width);
heightmatrix=maxrow-minrow;
widthmatrix=maxcol-mincol;
templatecols=ones(regionlength,1)*[0:height:(regionlength-1)*height];
templaterows=[1:regionlength]'*ones(1,regionlength);
templateindices=templatecols+templaterows;
YFP_raw(1)=NaN;
RFP_raw(1)=NaN;
linearindices=ones(regionlength*regionlength,1);
linearindices=linearindices*ones(1,numcells);
firstindex=sub2ind([height width],minrow,mincol);
for i=1:numcells
    indices=templateindices(1:heightmatrix(i),1:widthmatrix(i))+firstindex(i)-1;
    linearindices(1:numel(indices),i)=indices(:);
end
region_YFP_raw=YFP_raw(linearindices);
region_RFP_raw=RFP_raw(linearindices);
YFP_bg=prctile(region_YFP_raw,YFPprc)'; %adjust for cytoplasmic signal
RFP_bg=prctile(region_RFP_raw,RFPprc)';