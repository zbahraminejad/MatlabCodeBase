function DAs_pad=getnucmask(DAs_bs,nucr)

[DAs_ma,th1,th2,bg]=ThreshImage_MC_test(DAs_bs,0,25);  %10 bins higher for stricter threshold was log:60 or abs:30

DAs_ma=imfill(DAs_ma,'holes');
DAs_pad=padarray(DAs_ma,[1 1]);
DAs_bs=padarray(DAs_bs,[1 1]);
[height,width]=size(DAs_pad);
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));

%{
DAs_ma=DAs_pad(2:end-1,2:end-1);
DAs_bs(DAs_ma==0)=0;
cytogradients=edge(DAs_bs,'canny',[0.0001 0.003],1); %[.0063 .0156]
cytogradients=bwmorph(cytogradients,'thin',Inf);
cytogradients=imclose(cytogradients,strel('disk',2,0));
%cytogradients=imfill(cytogradients,'holes');
%cytogradients=imopen(cytogradients,strel('disk',1,0));
%}

%DAs_pad=padarray(DAs_ma,[1 1]);
%DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));  %this would remove stuff smaller than 25% of normal cell size.  Also removes spurs at edges & protruding points (but not intruding points, I think).
%fig1=figure; fig2=figure;
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
    [prnum,pcnum]=size(minimask);
    miniimg=DAs_bs(prmin:prmax,pcmin:pcmax);
    pixidx=find(minimask==1);
    pixints=miniimg(pixidx);
    pimin=min(pixints); pimax=max(pixints);
    kstep=0.02;
    [ksn,ksx]=ksdensity(pixints,pimin:kstep:pimax,'width',0.05);
    [wsn,wsx]=ksdensity(pixints,pimin:kstep:pimax,'width',0.1);
    downslope=find(diff(wsn)<0);
    firstpeak=downslope(1);
    upslope=find(diff(wsn)>0);
    upslope(upslope<=firstpeak)=[];
    if isempty(upslope)
        continue
    end
    nextvalley=upslope(1);
    downslope(downslope<=nextvalley)=[];
    nextpeak=downslope(1);
    downslope=find(diff(ksn)<0);
    kstep=(max(pixints)-min(pixints))/100;
    kbuffer=round(0.05/kstep);
    downslope(downslope>=nextpeak-kbuffer)=[];
    bgthresh=downslope(end);
    %}
    pixtodelete=pixidx(miniimg(pixidx)<ksx(bgthresh));
    if numel(pixtodelete)>=numel(pixidx)/2
        continue
    end
    
    %%% if objects separated, re-do bgsub on each object
    minimask(pixtodelete)=0;
    splitlabels=bwlabel(minimask);
    numsplitobs=max(max(splitlabels));
    if numsplitobs>1
        for mmi=1:numsplitobs
            splitpixidx=find(splitlabels==mmi);
            splitpixints=miniimg(splitpixidx);
            [ksn,ksx]=ksdensity(splitpixints,'width',0.05);
            firstpeak=find(diff(ksn)<0,1,'first');
            nextvalley=find(diff(ksn)>0);
            nextvalley(nextvalley<=firstpeak)=[];
            if isempty(nextvalley)
                continue
            end
            nextvalley=nextvalley(1);
            extrapixtodelete=splitpixidx(miniimg(splitpixidx)<ksx(nextvalley));
            if numel(extrapixtodelete)<numel(splitpixidx)/2
                pixtodelete=[pixtodelete;extrapixtodelete];
            end
        end
    end
    
    %%% remove pixtodelete from original mask
    [pixtodeleterow,pixtodeletecol]=ind2sub([prnum pcnum],pixtodelete);
    pixtodeleterow=pixtodeleterow+(prmin-1);
    pixtodeletecol=pixtodeletecol+(pcmin-1);
    pixtodelete=sub2ind([height width],pixtodeleterow,pixtodeletecol);
    DAs_pad(pixtodelete)=0;
end

DAs_pad=imfill(DAs_pad,'holes');
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));
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