function DAs_pad=getnucmask(DAs_bs,nucr)
clustersize=3000;
bgratio=1.5;
DAs_ma=ThreshImage_MC(DAs_bs,0);  %10 bins higher for stricter threshold was log:60 or abs:30

DAs_ma=imfill(DAs_ma,'holes');
DAs_pad=padarray(DAs_ma,[1 1]);
DAs_bs=padarray(DAs_bs,[1 1]);
[height,width]=size(DAs_pad);
DAs_pad=imopen(DAs_pad,strel('disk',round(nucr/2),0));

[nuclabels,obnum]=bwlabel(DAs_pad);
for ci=1:obnum       %ci = 4
    minimask=zeros(size(DAs_pad));
    pixidx=find(nuclabels==ci);
    minimask(pixidx)=1;
    %[pixr,pixc]=find(nuclabels==ci);
    [pixr,pixc]=ind2sub([height width],pixidx);
    prmin=min(pixr); prmax=max(pixr); pcmin=min(pixc); pcmax=max(pixc);
    %prnum=prmax-prmin+1; pcnum=pcmax-pcmin+1;
    minimask=minimask(prmin:prmax,pcmin:pcmax);
    minimaskorg=minimask;
    [prnum,pcnum]=size(minimask);
    tmap=zeros(prnum,pcnum);
    miniimg=DAs_bs(prmin:prmax,pcmin:pcmax);
    pixidx=find(minimask==1);
    pixints=miniimg(pixidx);
    bgthresh=findthreshold(pixints);
    if bgthresh==-10000
        continue
    end
    pixtodelete=pixidx(miniimg(pixidx)<bgthresh);
    if numel(pixtodelete)>=numel(pixidx)/bgratio
        continue
    end
    
    %%% if objects separated, re-do bgsub on each object %%%%%%%%%%%%%%%%%%
    minimask(pixtodelete)=0;
    minimask=imopen(minimask,strel('disk',round(nucr/2),0));
    [splitlabels,numsplitobs]=bwlabel(minimask);
    %numsplitobs=max(max(splitlabels));
    pixidx=find(splitlabels==1);
    if numsplitobs==1 && numel(pixidx)>clustersize
        pixints=miniimg(pixidx);
        bgthresh=findthreshold(pixints);
        if bgthresh==-10000
            continue
        end
        pixtoseparate=pixidx(miniimg(pixidx)<bgthresh);
        minimask(pixtoseparate)=0;
        tmap(pixtoseparate)=1;
        tmap=imopen(tmap,strel('disk',round(nucr*0.8),0));
        tmap=imerode(tmap,strel('disk',1,0));  %necessary to separate from remaining object(s)
        minimask=imopen(minimask,strel('disk',round(nucr/2),0));
        minimask(tmap==1)=1;
        [splitlabels,numsplitobs]=bwlabel(minimask);
    end
    if numsplitobs>1
        pixtodelete=[];
        for mmi=1:numsplitobs
            splitpixidx=find(splitlabels==mmi);
            splitpixints=miniimg(splitpixidx);
            bgthresh=findthreshold(splitpixints);
            if bgthresh==-10000
                continue
            end
            extrapixtodelete=splitpixidx(miniimg(splitpixidx)<bgthresh);
            if numel(extrapixtodelete)<numel(splitpixidx)/bgratio
                pixtodelete=[pixtodelete;extrapixtodelete];
            end
        end
        if numel(pixtodelete)>0
            minimask(pixtodelete)=0;
        end
    end
    deletedimage=minimaskorg-minimask;
    pixtodelete=find(deletedimage);
    
    %%% remove pixtodelete from original mask
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
    %{
    if ci==16   %32
        fprintf('ci = %0.0f\n',ci);
        pixints=DAs_bs(find(nuclabels==ci));
        pixidx=find(nuclabels==ci);
        tmap=zeros(size(DAs_pad));
        tmap(pixidx)=1;
        imshow(tmap);
        
        pixselect=pixidx(DAs_bs(pixidx)>th2);
        tmap=zeros(size(DAs_pad));
        tmap(pixselect)=1;
        imshow(tmap);
    end
    %}
    order=orderperimeter_nuc(coorset);
    orderedset=coorset(order,:);
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

function bgthresh = findthreshold(pixints)
ksthresh=-10000;
bgthresh=-10000;
pimin=min(pixints); pimax=max(pixints);
kstep=0.02;
[ksn,ksx]=ksdensity(pixints,pimin:kstep:pimax,'width',0.05);
upslope=diff(ksn)>0;
downslope=~upslope;
[downsegs,segnum]=bwlabel(downslope);
firstpeak=0;
for iseg=1:segnum
    segidx=find(downsegs==iseg);
    segidxmax=max(segidx);
    segintdiff=ksn(min(segidx))-ksn(segidxmax);
    segintdiffratio=segintdiff/ksn(segidxmax);
    if segintdiffratio>0.2  %real peak
        firstpeak=segidx(1);
        break
    end
end
[upsegs,segnum]=bwlabel(upslope);
for iseg=2:segnum
    segidx=find(upsegs==iseg);
    segidxmin=min(segidx);
    if segidxmin<firstpeak
        continue
    end
    segintdiff=ksn(max(segidx))-ksn(segidxmin);
    segintdiffper=segintdiff/ksn(segidxmin);
    if segintdiffper>0.2  %real 
        ksthresh=segidx(1);
        break
    end
end
if ksthresh~=-10000
    bgthresh=ksx(ksthresh);
end
end