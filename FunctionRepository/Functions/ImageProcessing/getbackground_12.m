function [primary_bg,secondary_bg]=getbackground_12(nuc_mask,centers,primary_raw,secondary_raw,nucr,winratio,resizeratio,primaryprc,secondaryprc)
nuc_mask=imresize(nuc_mask,resizeratio,'nearest');
primary_raw=imresize(primary_raw,resizeratio,'nearest');
secondary_raw=imresize(secondary_raw,resizeratio,'nearest');

offset=resizeratio*winratio*nucr;
regionlength=2*offset+1;
[height,width]=size(nuc_mask);
primary_raw(nuc_mask>0)=NaN;
secondary_raw(nuc_mask>0)=NaN;

centers=centers(~isnan(centers(:,1)),:);
centers=centers*resizeratio;
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
primary_raw(1)=NaN;
secondary_raw(1)=NaN;
linearindices=ones(regionlength*regionlength,1);
linearindices=linearindices*ones(1,numcells);
firstindex=sub2ind([height width],minrow,mincol);
for i=1:numcells
    indices=templateindices(1:heightmatrix(i),1:widthmatrix(i))+firstindex(i)-1;
    linearindices(1:numel(indices),i)=indices(:);
end
region_primary_raw=primary_raw(linearindices);
region_secondary_raw=secondary_raw(linearindices);
primary_bg=prctile(region_primary_raw,primaryprc)'; %adjust for cytoplasmic signal
secondary_bg=prctile(region_secondary_raw,secondaryprc)'; %adjust for cytoplasmic signal