function mask=product_getnucmask_rev03(DAs_bs,nucr)
%%% Perform nuclear segmentation by 'Concave Detection' & 'Histogram Segmentation' %%%%%%
% Author: Mingyu Chung
% Last Revision: 2/11/2013
%
% How to use this script:
% Arg1: background-subtracted image of cell nuclei (e.g. Hoechst stain)
% Arg2: average nuclear radius in pixels. This determines the minimum size
%       of nuclei to not consider as debris, and it only slightly matters
%       for segmenting the nuclei.
%
% What this script does:
% 1. Generates a mask purely based on intensity threshold
% 2. Segments nuclei by intensity distribution
% 3. Cleans up the mask by breaking some strands, spurs, etc.
% 4. For each object, searches perimeter and detects concave deflections
% 5. Bridges deflections/vertices to form segmentation borders
%    *Criteria applied for multiple vertices (e.g. more than 2 cells
%     overlapping, etc.)
% 
% This script requires:
% - segmentnuclei.m
% - bridge.m
%
% Users:
% Gautam Dey (gdey@stanford.edu)
% Steve Cappell (scappell@stanford.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucthird=round(nucr/3);

%%% define initial mask by threshold intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_bs=single(DAs_bs);
mask=getthreshold(DAs_bs,0);                     %2nd arg = # bins above threshold
mask=imfill(mask,'holes');

%%% segment mask by intensity distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=padarray(mask,[1 1]);                       %necessary for imopen
DAs_bs=padarray(DAs_bs,[1 1]);
mask=imopen(mask,strel('disk',nucthird,0)); %this would remove debris & spurs
mask=histogramsweep(mask,DAs_bs,nucthird);           %segment nuclei by intensity distribution

%%% clean up mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=imfill(mask,'holes');
mask=imopen(mask,strel('disk',nucthird,0));
mask=~bwmorph(~mask,'diag');                     %break diagonal connections
mask=~bwmorph(~mask,'bridge');                   %break orthogonal connections
nucedges=bwmorph(mask,'remove');                 %remove interior pixels
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(mask));

%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];                               %adjust to x-y convention
    [order,status]=getperimeter(coorset);
    if status==0                                 %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    bordermask=getborders(orderedset,bordermask,nucr);
end

%%% Incorporate segment borders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask(2:end-1,2:end-1);
bordermask=bordermask(2:end-1,2:end-1);
mask=mask & ~bordermask;                         %add segmentations
mask=~bwmorph(~mask,'diag');                     %remove diagonal connections
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ImageMask,th,bg]=getthreshold(orgimage,threshmargin)
TempSeries=orgimage(:);

%%% get background value: search distribution for most convex point %%%%%%%
tmax=max(TempSeries);
tmin=min(TempSeries);
nbin=200;tbin=(tmax-tmin)/nbin;
tmin=tmin+tbin/2;tmax=tmax-tbin/2;
[n,xout]=ksdensity(TempSeries,tmin:tbin:tmax);
gp=max([2,ceil(nbin/50)]);                      %gp=4
ng=getcurve(n,gp);                              %returns the angle change of intensity histogram over one bin at each point
[~,Ibg00]=min(ng);                              %find most convex point of intensity histogram
bg0=xout(Ibg00);

TempSeries2=TempSeries((TempSeries>(bg0-5*gp*tbin))&(TempSeries<(bg0+5*gp*tbin)));  %narrow search for peak
[n_bg,xout_bg]=ksdensity(TempSeries2);
[~,Ibg]=max(n_bg);                              %index of peak
bg = xout_bg(Ibg);

%%% get threshold: search upstream of background for most concave point %%%
upperbound = xout(Ibg00)+10*gp*tbin;
lowerbound = bg;
[n_th,xout_th]=ksdensity(TempSeries,lowerbound:(upperbound-bg)/100:upperbound);
ng_th=getcurve(n_th,gp);
pospeak=regionprops(bwlabel(ng_th>0),'PixelIdxList','Centroid','Area');    %same as negpeak, but for positive (concave) curvature
Ith00=zeros(1,size(pospeak,1));SC=Ith00;
for cc=1:size(pospeak,1)
    SC(cc)=sum(ng_th(pospeak(cc).PixelIdxList));                           %sum of curvatures for each pospeak group
    Ith00(cc)=round(pospeak(cc).Centroid(1));                              %index of center of each pospeak group
end
[~,I]=max(SC);
th=xout_th(Ith00(I)+threshmargin);                                         %I added 10 so the threshold is higher

%%% get mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImageMask=single(orgimage>th);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ng=getcurve(n,gp)
ns=[n(1+gp:-1:1+1),n,n(end-1:-1:end-gp)];
theta1=atan2(ns(1+gp:end-gp)-ns(1:end-2*gp),gp);
theta2=atan2(ns(1+2*gp:end)-ns(1+gp:end-gp),gp);
ng_a=acos(cos(theta2-theta1));
ng_v=ns(1:end-2*gp)+ns(1+2*gp:end)-2*ns(1+gp:end-gp);
ng=ng_a.*(ng_v>0)-ng_a.*(ng_v<0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DAs_pad=histogramsweep(DAs_pad,DAs_bs,nucthird)
clustersize=500;
bgratio=1.5;
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
    bgthresh=findbestthreshold(pixints);
    if bgthresh==-10000
        continue
    end
    pixtodelete=pixidx(miniimg(pixidx)<bgthresh);
    minimask(pixtodelete)=0;
    minimask=imopen(minimask,strel('disk',nucthird,0));
    
    %%% repeat process until no further segmentation occurs %%%%%%%%%%%%%%%
    [splitlabels,numsplitobs]=bwlabel(minimask);
    mmi=1;
    while mmi<=numsplitobs
        pixidx=find(splitlabels==mmi);
        if numel(pixidx)>clustersize
            separatemap=zeros(prnum,pcnum);
            remainmap=separatemap;
            pixints=miniimg(pixidx);
            bgthresh=findbestthreshold(pixints);
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
            separatemap=imopen(separatemap,strel('disk',nucthird,0));
            separatemap=imerode(separatemap,strel('disk',1,0));  %necessary to separate from remaining object(s)
            pixtokeep=pixidx(~idxtoseparate);
            remainmap(pixtokeep)=1;
            remainmap=imopen(remainmap,strel('disk',nucthird,0));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bgthresh = findbestthreshold(pixints)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [order,status]=getperimeter(set)
numpoints=size(set,1);
status=1; %change to zero if bad
order=zeros(numpoints,1);
order(1)=1;
%%% determine clockwise direction from r(1),c(1) %%%%%%%%%%%%%%%%%%
compass=[1,0;1,-1;0,-1;-1,-1;-1,0;-1,1;0,1;1,1];
buffer=zeros(size(set));
x=set(1,1); y=set(1,2);
buffer(:,1)=set(:,1)-x;
buffer(:,2)=set(:,2)-y;
buffer(1,:)=[10000,10000];
adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
[~,idx]=sortrows(buffer(adjidx,:),[-2 1]);   %sort first by descending y then ascending x
adjidx=adjidx(idx(1),:);
%%% order ring coordinates contiguously %%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=2:numpoints-1
    order(pp)=adjidx;
    x=buffer(adjidx,1); y=buffer(adjidx,2);
    buffer(:,1)=buffer(:,1)-x;
    buffer(:,2)=buffer(:,2)-y;
    buffer(adjidx,:)= [10000,10000];
    adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
    if numel(adjidx)>1
        pole=find(compass(:,1)==x & compass(:,2)==y);
        antipole=pole+4;
        priority=[antipole+1:antipole+7]';
        priority=mod(priority,8); priority(priority==0)=8;
        priority=compass(priority,:);
        poles=find(ismember(priority,buffer(adjidx,:),'rows'));
        x=priority(poles(1),1); y=priority(poles(1),2);
        adjidx=find(buffer(:,1)==x & buffer(:,2)==y);
    elseif numel(adjidx)==0                 %error: kickout with status flag
        status=0;
        break
    end
end
if status==1
    order(numpoints)=adjidx;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bordermask=getborders(orderedset,bordermask,nucr)
bridgeimage=0;
nucr=round(nucr/4)*4;
%minperi=50;  %minimum length of perimeter to consider breaking

%%% calculate tangent angle & detect vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
offsetshort=nucr/4;
offsetlong=2*offsetshort;
gradientoffsetshort=offsetshort;
gradientoffsetlong=offsetlong;

orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
shortdiff=orderedsetoffsetshort-orderedset;
longdiff=orderedsetoffsetlong-orderedset;
shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians
longgrad=atan2(longdiff(:,2),longdiff(:,1));
shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
shortgraddiff=shortgradoffset-shortgrad;
longgraddiff=longgradoffset-longgrad;
shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);
longgraddiff=longgraddiff+2*pi*(longgraddiff<0);
shortgradthresh = pi/6;
longgradthresh = pi/6;
vIdxmasklong=longgraddiff>longgradthresh & longgraddiff<pi;
vIdxmasklong=[zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
vIdxmasklong=imdilate(vIdxmasklong,strel('square',1+nucr/2));
shortgraddiff(shortgraddiff>=pi)=0;
vIdxmaskshort=shortgraddiff>shortgradthresh;

vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));
vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
maxmask=zeros(size(vIdxmaskshort));
for rpc=1:length(vIdxobs)
    pix=vIdxobs(rpc).PixelIdxList;
    [~,index]=max(shortgraddiff(pix));
    maxmask(pix(index)+offsetshort)=1;
end

vIdxmask=vIdxmasklong & maxmask;
vIdx=find(vIdxmask);

%%% mark candidate vertices [debugging purposes] %%%%%%%%%%%%%%%%%%%%%%%%%%
if bridgeimage
    vpos=orderedset(vIdx,:);
    for vc=1:size(vpos,1)
        bordermask(vpos(vc,2),vpos(vc,1))=1;
    end
end

%%% pair and connect vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vnum=length(vIdx);
if vnum>=2
    periIdx=vIdx;
    perisize=length(vIdxmask);
    periIdxadj1=[periIdx(2:end);perisize+periIdx(1)];
    pairperi1=periIdxadj1-periIdx;    %perimeter distance between adj vertices
end
while vnum>=2
    skipvertices=0;
    vpos=orderedset(vIdx,:);
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    [bestpair1,ordx1]=sort(pairperi1./pairdist1);
    idx=ordx1(end);
    if idx==vnum
        idxadj=1;
    else
        idxadj=idx+1;
    end
    if vnum>=5
        vposadj2=[vposadj1(2:end,:);vposadj1(1,:)];
        pair2=vposadj2-vpos;
        pairdist2=sqrt(sum(pair2.^2,2));
        pairperi2=[pairperi1(2:end);pairperi1(1)];
        pairperi2=pairperi1+pairperi2;    %perimeter btwn every other vertice
        [bestpair2,ordx2]=sort(pairperi2./pairdist2);
        if bestpair2(end)>bestpair1(end)
            skipvertices=1;
            idx=ordx2(end);
            if idx==vnum
                idxadj=2; idxint=1;
            elseif idx==vnum-1
                idxadj=1; idxint=vnum;
            else
                idxadj=idx+2; idxint=idx+1;
            end
        end
    end

    if vnum<=3
        if vnum==3
            idxadjadj=idxadj+1;
            if idxadjadj>vnum
                idxadjadj=1;
            end
            pairperi1(idxadj)=pairperi1(idxadj)+pairperi1(idxadjadj);
            pairperi1(idxadjadj)=nucr/2*pi;
        end
        vlo=pairperi1<nucr/2*pi;   %less than half circumference of avg size nucleus
        if sum(vlo)             %either side too short
            break
        end
    end

    [bcx,bcy]=buildbridge(vpos(idx,:),vpos(idxadj,:));
    for bci=1:length(bcx)
        bordermask(bcy(bci),bcx(bci))=1;
    end
    
    %%% assign new perimter distance & remove old vertices %%%%%%%%%%%%%%%%
    previdx=idx-1;
    if previdx==0
        previdx=vnum;
    end
    pairperi1(previdx)=pairperi1(previdx)+length(bcx)-1+pairperi1(idxadj);
    if skipvertices==0
        vIdx([idx,idxadj])=[];
        pairperi1([idx,idxadj])=[];
    else
        vIdx([idx,idxint,idxadj])=[];
        pairperi1([idx,idxint,idxadj])=[];
    end
    vnum=length(vIdx);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bx,by] = buildbridge(vpos1,vpos2)
lengthx=vpos2(1)-vpos1(1);
lengthy=vpos2(2)-vpos1(2);
longerside=max([abs(lengthx) abs(lengthy)]);
stepx=lengthx/longerside;
stepy=lengthy/longerside;
bx=zeros(longerside+1,1);by=zeros(longerside+1,1);
for bs=0:longerside
    bx(1+bs)=vpos1(1)+round(bs*stepx);
    by(1+bs)=vpos1(2)+round(bs*stepy);
end
end