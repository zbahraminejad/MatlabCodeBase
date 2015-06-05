function mask=segmentdeflections_singlecell(mask,nucr,cellsizeratio)
%%% only consider big nuclei %%%%%
bigthresh=round(cellsizeratio*pi*nucr^2);  %1.33*pi*nucr^2 is the avg cell size
onlybigs=bwareaopen(mask,bigthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucedges=bwmorph(onlybigs,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(mask));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];  %adjust to x-y convention
    [order,status]=orderperimeter([c,r]);
    if status==0    %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    bordermask=splitdeflections(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,round(0.33*pi*nucr^2));
end