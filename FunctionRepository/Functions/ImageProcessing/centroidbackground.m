function nuc_bg=centroidbackground(nuc_mask,centers,raw,nucr,winratio,resizeratio,prcthresh)
nuc_mask=imresize(nuc_mask,resizeratio,'nearest');
raw=imresize(raw,resizeratio,'nearest');

offset=resizeratio*winratio*nucr;
regionlength=2*offset+1;
[height,width]=size(nuc_mask);
raw(nuc_mask>0)=NaN;

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
raw(1)=NaN;
linearindices=ones(regionlength*regionlength,1);
linearindices=linearindices*ones(1,numcells);
firstindex=sub2ind([height width],minrow,mincol);
for i=1:numcells
    indices=templateindices(1:heightmatrix(i),1:widthmatrix(i))+firstindex(i)-1;
    linearindices(1:numel(indices),i)=indices(:);
end
region_raw=raw(linearindices);
nuc_bg=prctile(region_raw,prcthresh)';

%%% debug %%%%%%%%%
%{
threshprctile=ones(numcells,1)*NaN;
for i=1:size(region_raw,2)
    threshprctile(i)=firstpeak(region_raw(:,i));
end
hist(threshprctile,100);

values=region_raw(:,i);
kmin=min(values); kmax=max(values); kstep=(kmax-kmin)/100;
[ksn,ksx]=ksdensity(values,kmin:kstep:kmax);
ksthresh=find(diff(ksn)<0,1,'first');
ksthresh=(ksthresh-1);
plot(ksx,ksn)
fprintf('thresh percentile = %0.0f\n',ksthresh);
fprintf('10th percentile = %0.0f\n',prctile(values,10));
fprintf('25th percentile = %0.0f\n',prctile(values,25));
fprintf('30th percentile = %0.0f\n',prctile(values,30));
fprintf('50th percentile = %0.0f\n',prctile(values,50));
%}
