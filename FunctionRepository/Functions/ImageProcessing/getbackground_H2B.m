function nuc_bg=getbackground_H2B(nuc_mask,centers,nuc_raw,nucr,ratio)
nuc_mask=imresize(nuc_mask,ratio,'nearest');
nuc_raw=imresize(nuc_raw,ratio,'nearest');

offset=ratio*15*nucr;
regionlength=2*offset+1;
[height,width]=size(nuc_mask);
nuc_raw(nuc_mask>0)=NaN;

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
nuc_raw(1)=NaN;
linearindices=ones(regionlength*regionlength,1);
linearindices=linearindices*ones(1,numcells);
firstindex=sub2ind([height width],minrow,mincol);
for i=1:numcells
    indices=templateindices(1:heightmatrix(i),1:widthmatrix(i))+firstindex(i)-1;
    linearindices(1:numel(indices),i)=indices(:);
end
region_nuc_raw=nuc_raw(linearindices);
nuc_bg=nanmedian(region_nuc_raw)';
