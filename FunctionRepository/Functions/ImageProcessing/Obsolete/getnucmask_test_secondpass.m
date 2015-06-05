function DAs_pad=getnucmask(DAs_bs,nucr)

DAs_ma=ThreshImage_MC(DAs_bs,0);  %10 bins higher for stricter threshold
DAs_ma=imfill(DAs_ma,'holes');

%DAs_pad=zeros(size(DAs_ma,1)+2,size(DAs_ma,2)+2);
%DAs_pad(2:end-1,2:end-1)=DAs_ma;
DAs_pad=padarray(DAs_ma,[1 1]);
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));  %this would remove stuff smaller than 25% of normal cell size.  Also removes spurs at edges & protruding points (but not intruding points, I think).

%DAs_pad=imdilate(DAs_pad,strel('disk',2,0));    %this will actually be used for cytoline
%DAs_pad=imfill(DAs_pad,'holes');
%DAs_pad=~bwmorph(~DAs_pad,'diag');         %break connections
%DAs_pad=~bwmorph(~DAs_pad,'bridge');      %break connections
[DAs_la,obnum]=bwlabel(DAs_pad);
DAs_rp=regionprops(DAs_la,'Area','PixelIdxList');
%ringlabels=DAs_la-imerode(DAs_la,strel('disk',1,0));

%%% refine borders for big blobs %%%%%%%%%%%%%%%%%%%%%%
[~,nextthresh,~]=ThreshImage_MC(DAs_bs,10);
bigthresh=nucr*nucr*3*6;
for ci=1:obnum
    if DAs_rp(ci).Area>bigthresh
        pixidx=DAs_rp(ci).PixelIdxList;
        pixvals=DAs_bs(pixidx);
        %lowthresh=prctile(pixvals,25)-1.5*iqr(pixvals);
        lowpix=pixidx(pixvals<nextthresh);
        DAs_pad(lowpix)=0;
        %{
        [lowthresh,isbg]=ThreshImage_MC_bgcheck(pixvals,superbg,0);
        if isbg
            lowpix=pixidx(pixvals<lowthresh);
            DAs_pad(lowpix)=0;
        end
        %}
    end
end

%%% clarify borders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));
DAs_pad=imfill(DAs_pad,'holes');
DAs_pad=~bwmorph(~DAs_pad,'diag');         %break connections
DAs_pad=~bwmorph(~DAs_pad,'bridge');      %break connections
nucedges=bwmorph(DAs_pad,'remove');
[ringlabels,obnum]=bwlabel(nucedges);

%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bordermask=zeros(size(DAs_pad));
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