function DAs_pad=histsweep(DAs_pad,DAs_bs,nucfifth)
clustersize=500;
bgratio=1.5; %1.5
[height,width]=size(DAs_pad);
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
    %%%%%%% avoid removing nucleoli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minimask_holes=minimask;
    minimask_holes(pixtodelete)=0;
    minimask=imfill(minimask_holes,'holes');
    minimask_holes= minimask_holes==0 & minimask==1;
    pixholes=find(minimask_holes==1);
    miniimg(pixholes)=1.5;  %this is the threshold intensity above which things will not be removed
    minimask=imopen(minimask,strel('disk',nucfifth,0));
    
    %%% repeat process until no further segmentation occurs %%%%%%%%%%%%%%%
    [splitlabels,numsplitobs]=bwlabel(minimask);
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
            if numel(pixtoseparate)>numel(pixidx)/bgratio   %don't remove because it's larger than the 'foreground'
                mmi=mmi+1;
                continue
            end
            separatemap(pixtoseparate)=1;
            separatemap=imopen(separatemap,strel('disk',nucfifth,0));  %ignore cells that are too small
            separatemap=imerode(separatemap,strel('disk',1,0));  %necessary to separate from remaining object(s)
            pixtokeep=pixidx(~idxtoseparate);
            remainmap(pixtokeep)=1;
            %%%%%%%% avoid removing nucleoli %%%%%%%%%%%%%%%%%%%%%%%
            remainmap_holes=remainmap;
            remainmap=imfill(remainmap_holes,'holes');
            remainmap_holes=remainmap_holes==0 & remainmap==1;
            pixholes=find(remainmap_holes==1);
            miniimg(pixholes)=1.5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            remainmap=imopen(remainmap,strel('disk',nucfifth,0));
            [rmlabel,rmnum]=bwlabel(remainmap);
            %if ismember(1,separatemap) || rmnum>1  %segmentation occurred
                rmlabel=rmlabel+numsplitobs;
                rmlabel(rmlabel==numsplitobs)=0;
                numsplitobs=numsplitobs+rmnum;
            %end
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
    segidxmin=segidx(1);
    segidxmax=segidx(end);
    segintdiff=ksn(segidxmin)-ksn(segidxmax);
    segintdiffratio=segintdiff/ksn(segidxmax);
    if segintdiffratio>0.1  %real peak
        ksthresh=segidxmax+1;
        break
    end
end
if ksthresh<midstep
    bgthresh=ksx(ksthresh);
end
end