function DAs_pad=getnucmask(DAs_bs,nucr,bigratio)
debrisrad=round(nucr/6);
debrisarea=round(pi*(nucr/4)^2);
DAs_ma=ThreshImage_MC(DAs_bs,0);
DAs_ma=imfill(DAs_ma,'holes');
DAs_pad=padarray(DAs_ma,[1 1]);
DAs_pad=imopen(DAs_pad,strel('disk',debrisrad,0));
DAs_pad=bwareaopen(DAs_pad,debrisarea);      %added this to do stronger debris sweep without removing elliptical nuclei
DAs_pad=~bwmorph(~DAs_pad,'diag');      %break connections
DAs_pad=~bwmorph(~DAs_pad,'bridge');    %break connections
%%% only consider big nuclei %%%%%
bigthresh=round(bigratio*pi*nucr^2);
onlybigs=bwareaopen(DAs_pad,bigthresh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucedges=bwmorph(onlybigs,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(DAs_pad));
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
    bordermask=segmentnuclei_skipper(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad(2:end-1,2:end-1);
bordermask=bordermask(2:end-1,2:end-1);
DAs_pad=DAs_pad & ~bordermask;
DAs_pad=~bwmorph(~DAs_pad,'diag');
end