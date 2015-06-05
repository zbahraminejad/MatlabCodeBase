function DAs_pad=getnucmask(DAs_bs,nucr)

DAs_ma=ThreshImage_MC(DAs_bs,0);  %10 bins higher for stricter threshold
DAs_ma=imfill(DAs_ma,'holes');

DAs_pad=padarray(DAs_ma,[1 1]);
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));  %this would remove stuff smaller than 25% of normal cell size.  Also removes spurs at edges & protruding points (but not intruding points, I think).
DAs_pad=~bwmorph(~DAs_pad,'diag');         %break connections
DAs_pad=~bwmorph(~DAs_pad,'bridge');      %break connections
nucedges=bwmorph(DAs_pad,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(DAs_pad));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    set=[c,r];  %adjust to x-y convention
    order=orderperimeter_nuc(set);
    orderedset=set(order,:);
    bordermask=segmentnuclei(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad(2:end-1,2:end-1);
%bordermask=imdilate(bordermask,strel('disk',1,8)); %required because bwlabel can't discriminate a one-pixel separation
%bordermask=bwmorph(bordermask,'diag');
bordermask=bordermask(2:end-1,2:end-1);
DAs_pad=DAs_pad & ~bordermask;
DAs_pad=~bwmorph(~DAs_pad,'diag');
%DAs_pad=~bwmorph(~DAs_pad,'diag'); %clear remaining strands
end