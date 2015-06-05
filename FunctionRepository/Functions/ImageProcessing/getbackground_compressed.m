function [nuc_bs,YFP_bs,RFP_bs]=getbackground(nuc_mask,nuc_info,nuc_raw,YFP_raw,RFP_raw,nucr,ratio)
nuc_mask=imresize(nuc_mask,ratio,'nearest');
nuc_raw=imresize(nuc_raw,ratio,'nearest');
YFP_raw=imresize(YFP_raw,ratio,'nearest');
RFP_raw=imresize(RFP_raw,ratio,'nearest');

offset=ratio*10*nucr;
numcells=length(nuc_info);
[height,width]=size(nuc_mask);
nuc_raw(nuc_mask>0)=NaN;
YFP_raw(nuc_mask>0)=NaN;
RFP_raw(nuc_mask>0)=NaN;
nuc_bs=zeros(numcells,1);
YFP_bs=zeros(numcells,1);
RFP_bs=zeros(numcells,1);
for i=1:numcells
    center=round(nuc_info(i).Centroid*ratio);
    minrow=center(2)-offset; minrow=max([minrow 1]);
    maxrow=center(2)+offset; maxrow=min([maxrow height]);
    mincol=center(1)-offset; mincol=max([mincol 1]);
    maxcol=center(1)+offset; maxcol=min([maxcol width]);
    nucwindow=nuc_raw(minrow:maxrow,mincol:maxcol);
    YFPwindow=YFP_raw(minrow:maxrow,mincol:maxcol);
    RFPwindow=RFP_raw(minrow:maxrow,mincol:maxcol);
    nuc_bs(i)=nanmedian(nucwindow(:));
    YFP_bs(i)=prctile(YFPwindow(:),10); %adjust for cytoplasmic signal
    RFP_bs(i)=nanmedian(RFPwindow(:));
    %imshow(imadjust(mat2gray(YFP_raw(minrow:maxrow,mincol:maxcol))));
end
