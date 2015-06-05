function DAs_pad=getnucmask_histoseg(DAs_bs,nucr)
clustersize=1000;
bgratio=1.5;
DAs_ma=ThreshImage_MC(DAs_bs,0);  %10 bins higher for stricter threshold was log:60 or abs:30

DAs_ma=imfill(DAs_ma,'holes');
DAs_pad=padarray(DAs_ma,[1 1]);
DAs_bs=padarray(DAs_bs,[1 1]);
[height,width]=size(DAs_pad);
DAs_pad=imopen(DAs_pad,strel('disk',round(nucr/2),0));

[nuclabels,obnum]=bwlabel(DAs_pad);
for ci=1:obnum
    minimask=zeros(size(DAs_pad));
    pixidx=find(nuclabels==ci);
    minimask(pixidx)=1;
    [pixr,pixc]=ind2sub([height width],pixidx);
    prmin=min(pixr); prmax=max(pixr); pcmin=min(pixc); pcmax=max(pixc);
    minimask=minimask(prmin:prmax,pcmin:pcmax);
    minimaskorg=minimask;
    [prnum,pcnum]=size(minimask);
    miniimg=DAs_bs(prmin:prmax,pcmin:pcmax);
    
    %%% remove background pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pixidx=find(minimask==1);
    pixints=miniimg(pixidx);
    bgthresh=findthreshold(pixints);
    if bgthresh==-10000
        continue
    end
    pixtodelete=pixidx(miniimg(pixidx)<bgthresh);
    minimask(pixtodelete)=0;
    minimask=imopen(minimask,strel('disk',round(nucr/2),0));
    
    %%% repeat process until no further segmentation occurs %%%%%%%%%%%%%%%
    [splitlabels,numsplitobs]=bwlabel(minimask);
    %for mmi=1:numsplitobs
    mmi=1;
    while mmi<=numsplitobs
        pixidx=find(splitlabels==mmi);
        if numel(pixidx)>clustersize
            separatemap=zeros(prnum,pcnum);
            remainmap=separatemap;
            pixints=miniimg(pixidx);
            bgthresh=findthreshold(pixints);
            if bgthresh==-10000 || bgthresh>1.5
                mmi=mmi+1;
                continue
            end

            idxtoseparate=miniimg(pixidx)<bgthresh;
            pixtoseparate=pixidx(idxtoseparate);
            if numel(pixtoseparate)>numel(pixidx)/bgratio
                mmi=mmi+1;
                continue
            end
            separatemap(pixtoseparate)=1;
            separatemap=imopen(separatemap,strel('disk',round(nucr/2),0));
            separatemap=imerode(separatemap,strel('disk',1,0));  %necessary to separate from remaining object(s)
            pixtokeep=pixidx(~idxtoseparate);
            remainmap(pixtokeep)=1;
            remainmap=imopen(remainmap,strel('disk',round(nucr/2),0));
            [rmlabel,rmnum]=bwlabel(remainmap);
            if ismember(1,separatemap) || rmnum>1  %segmentation occurred
                rmlabel=rmlabel+numsplitobs;
                rmlabel(rmlabel==numsplitobs)=0;
                numsplitobs=numsplitobs+rmnum;
            end
            splitlabels(pixidx)=0;
            splitlabels(separatemap==1)=1;
            splitlabels=splitlabels+rmlabel;
        end
        mmi=mmi+1;
    end
    minimask=splitlabels>0;
    deletedimage=minimaskorg-minimask;
    pixtodelete=find(deletedimage);
    
    %%% remove pixtodelete from original mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pixtodeleterow,pixtodeletecol]=ind2sub([prnum pcnum],pixtodelete);
    pixtodeleterow=pixtodeleterow+(prmin-1);
    pixtodeletecol=pixtodeletecol+(pcmin-1);
    pixtodelete=sub2ind([height width],pixtodeleterow,pixtodeletecol);
    DAs_pad(pixtodelete)=0;
end

DAs_pad=imfill(DAs_pad,'holes');
DAs_pad=imopen(DAs_pad,strel('disk',round(nucr/2),0));
DAs_pad=~bwmorph(~DAs_pad,'diag');         %break connections
DAs_pad=~bwmorph(~DAs_pad,'bridge');      %break connections
nucedges=bwmorph(DAs_pad,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(DAs_pad));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];  %adjust to x-y convention
    order=orderperimeter_nuc(coorset);
    orderedset=coorset(order,:);
    bordermask=segmentnuclei(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad(2:end-1,2:end-1);
bordermask=bordermask(2:end-1,2:end-1);
DAs_pad=DAs_pad & ~bordermask;
DAs_pad=~bwmorph(~DAs_pad,'diag');
end

function bgthresh = findthreshold(pixints)
ksthresh=10000;
bgthresh=-10000;
pimin=min(pixints); pimax=max(pixints);
kstep=0.02;
midstep=round((pimax-pimin)/(kstep*2));
[ksn,ksx]=ksdensity(pixints,pimin:kstep:pimax,'width',0.05);
downslope=diff(ksn)<0;
[downsegs,segnum]=bwlabel(downslope);
for iseg=1:segnum
    segidx=find(downsegs==iseg);
    segidxmax=max(segidx);
    segintdiff=ksn(min(segidx))-ksn(segidxmax);
    segintdiffratio=segintdiff/ksn(segidxmax);
    if segintdiffratio>0.2  %real peak
        ksthresh=segidx(end)+1;
        break
    end
end
if ksthresh<midstep
    bgthresh=ksx(ksthresh);
end
end