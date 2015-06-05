function [nuc_bs,YFP_bs,RFP_bs]=getbackground(nuc_mask,nuc_info,nuc_raw,YFP_raw,RFP_raw,nucr,ratio)
nuc_mask=imresize(nuc_mask,ratio,'nearest');
nuc_raw=imresize(nuc_raw,ratio,'nearest');
YFP_raw=imresize(YFP_raw,ratio,'nearest');
RFP_raw=imresize(RFP_raw,ratio,'nearest');

offset=ratio*10*nucr;
regionlength=2*offset+1;
numcells=length(nuc_info);
[height,width]=size(nuc_mask);
nuc_raw(nuc_mask>0)=NaN;
YFP_raw(nuc_mask>0)=NaN;
RFP_raw(nuc_mask>0)=NaN;

nuc_info=struct2cell(nuc_info');
center=squeeze(cell2mat(nuc_info(2,1,:)))*ratio;
x=center(1,:); y=center(2,:);

onesmatrix=ones(1,numcells);
minrow=round(y-offset); minrow=max(minrow,onesmatrix);
maxrow=round(y+offset); maxrow=min(maxrow,onesmatrix*height);
mincol=round(x-offset); mincol=max(mincol,onesmatrix);
maxcol=round(x+offset); maxcol=min(maxcol,onesmatrix*width);
heightmatrix=maxrow-minrow;
widthmatrix=maxcol-mincol;
templatecols=ones(regionlength,1)*[0:height:(regionlength-1)*height];
templaterows=[1:regionlength]'*ones(1,regionlength);
templateindices=templatecols+templaterows;
nuc_raw(1)=NaN; YFP_raw(1)=NaN; RFP_raw(1)=NaN; %will use index 1 as dummy NaN
linearindices=ones(regionlength*regionlength,1);
linearindices=linearindices*ones(1,numcells);
firstindex=sub2ind([height width],minrow,mincol);
for i=1:numcells
    indices=templateindices(1:heightmatrix(i),1:widthmatrix(i))+firstindex(i)-1;
    linearindices(1:numel(indices),i)=indices(:);
end
region_nuc_raw=nuc_raw(linearindices);
region_YFP_raw=YFP_raw(linearindices);
region_RFP_raw=RFP_raw(linearindices);
nuc_bs=nanmedian(region_nuc_raw);
YFP_bs=prctile(region_YFP_raw,10); %adjust for cytoplasmic signal
RFP_bs=nanmedian(region_RFP_raw);
