function mask=segmentdeflections(mask,nucr,cellsizeratio)
%%% only consider big nuclei %%%%%
bigthresh=round(cellsizeratio*pi*nucr^2);  %1.33*pi*nucr^2 is the avg cell size
onlybig=bwareaopen(mask,bigthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucedges=bwmorph(onlybig,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(mask));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    numpoints=length(r);
    if numpoints<20    %no need to segment
        continue
    end
    coorset=[c,r];                      %adjust to x-y convention
    [order,status]=orderperimeter(coorset);
    if status==0                    %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    bordermask=splitdeflections_test(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,round(0.33*pi*nucr^2)); %previously 0.33
end