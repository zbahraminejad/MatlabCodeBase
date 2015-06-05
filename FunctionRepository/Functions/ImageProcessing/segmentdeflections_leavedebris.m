function mask=segmentdeflections(mask,nucr)
nucedges=bwmorph(mask,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(mask));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];                      %adjust to x-y convention
    [order,status]=orderperimeter(coorset);
    if status==0                    %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    %bordermask=splitdeflections_3(orderedset,bordermask,nucr);
    bordermask=splitdeflections_4(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
end