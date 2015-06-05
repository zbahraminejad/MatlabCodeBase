function DAs_pad=getnucmaskandring(DAs_bs,nucr)

DAs_ma=ThreshImage_MC(DAs_bs,0);  %10 bins higher for stricter threshold
DAs_ma=imfill(DAs_ma,'holes');

DAs_pad=zeros(size(DAs_ma,1)+2,size(DAs_ma,2)+2);
DAs_pad(2:end-1,2:end-1)=DAs_ma;
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2));  %this would remove stuff smaller than 25% of normal cell size.  Also removes spurs at edges & protruding points (but not intruding points, I think).
DAs_pad=bwmorph(DAs_pad,'fill');
nucedges=bwmorph(DAs_pad,'diag');   %removes 8-connectivity of background (e.g. single pixel connections)
nucedges=bwmorph(nucedges,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(DAs_pad));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    %fprintf('ci = %0.0f\n',ci);
    [r,c]=find(ringlabels==ci);
    rlidx=find(ringlabels==ci);
    set=[c,r];  %adjust to x-y convention
    bordermask=segmentnucleiandring(set,rlidx,DAs_pad,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad(2:end-1,2:end-1);
%bordermask=imdilate(bordermask,strel('disk',1,8)); %required because bwlabel can't discriminate a one-pixel separation
bordermask=bwmorph(bordermask,'diag');
bordermask=bordermask(2:end-1,2:end-1);
DAs_pad=DAs_pad & ~bordermask;
DAs_pad=~bwmorph(~DAs_pad,'diag');
end