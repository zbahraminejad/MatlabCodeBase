function mask=segmentmerges(centers,centersprev,masses,massesprev,orgmask,nucr)
%%% detect merges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(orgmask);
winrad=8*nucr;
xcur=centers(:,1); ycur=centers(:,2);
xprev=centersprev(:,1); yprev=centersprev(:,2);
numcells=numel(xcur);
merge=zeros(numcells,1);
for i=1:numcells
    neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    massdiff=(masses(i)-massesprev(neighbors));
    [~,closestidx]=sort(massdiff);
    samecell=closestidx(1);
    merge(i)=massdiff(samecell)<0.2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucedges=bwmorph(onlybigs,'remove');
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
    bordermask=segmentnuclei_skipper(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
end