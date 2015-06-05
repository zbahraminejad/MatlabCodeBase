function newmask=segmentdeflections_parallel(mask,nucr,cellsizeratio)
bordermask=zeros(size(mask));
sizethresh=round(cellsizeratio*pi*nucr^2);  %1.33*pi*nucr^2 is the avg cell size
gatedmask=bwareaopen(mask,sizethresh);
nucedge=bwmorph(gatedmask,'remove');
%[ringlabels,obnum]=bwlabel(nucedge);
nucedgeidx=struct2cell(regionprops(nucedge,'PixelIdsList'));



    [r,c]=find(ringlabels==ci);
    numpoints=length(r);
    coorset=[c,r];                      %adjust to x-y convention
    [order,status]=orderperimeter(coorset);
    if status==0                    %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    
    bordermask=segmentnuclei_skipper(orderedset,bordermask,nucr);






mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
end