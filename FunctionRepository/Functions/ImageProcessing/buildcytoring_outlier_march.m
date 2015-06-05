function [DAs_da,finalcytoring]=buildcytoring_outlier_march(DAs_pad,REs_bs,nucr)


if nucr==16
    outerrad=4;
    midrad=3;
    innerrad=2;
elseif nucr==8
    outerrad=3;  %was 3
    midrad=2;
    innerrad=2;
end
zerorad=innerrad-1;
ringwidth=outerrad-zerorad;
%ringmargin=outerrad+2;
%}
[height,width]=size(DAs_pad);
zeroborder=ones(height,width);
zeroborder(1,:)=0;zeroborder(end,:)=0;zeroborder(:,1)=0;zeroborder(:,end)=0;
ringzeroborder=ones(height,width);
ringzeroborder(1:ringwidth+1,:)=0;ringzeroborder(end-ringwidth:end,:)=0;ringzeroborder(:,1:ringwidth+1)=0;ringzeroborder(:,end-ringwidth:end)=0;

borders=bwmorph(DAs_pad,'bothat');
borders=imdilate(borders,strel('disk',ringwidth,8));

%%%%%%%%%%%%%%%% new addition %%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad & zeroborder;   %necessary to imerode from edges
realnuc_la=bwlabel(DAs_pad);
smoothmask=imopen(DAs_pad,strel('disk',nucr/2,0));
realnuc_la=realnuc_la.*smoothmask;  %maintains IDs, but smoothens & removes debris
DAs_da=regionprops(realnuc_la,'Area','Centroid','PixelIdxList');
numcells=length(DAs_da);
%DAs_cen=regionprops(realnuc_la,'Centroid');
%centroidlist=cell2mat(struct2cell(DAs_cen))';
crm_la_org=bwlabel(bwmorph(realnuc_la,'thicken',midrad));
crm_la=zeros(height,width);
for i=1:numcells
    if DAs_da(i).Area==0
        continue
    end
    x=round(DAs_da(i).Centroid(1));
    y=round(DAs_da(i).Centroid(2));
    crm_label=crm_la_org(y,x);
    if crm_label>0
        crm_la(crm_la_org==crm_label)=i;
    end
end
smoothmask=imopen(crm_la,strel('disk',nucr/4,0));
crm_la=crm_la.*logical(smoothmask);  %maintains IDs, but smoothens & removes debris
cro_la=imdilate(crm_la,strel('disk',1,0));
cri_la=imerode(crm_la,strel('disk',1,0));
cytoline=crm_la-cri_la;
cytolineskel=bwmorph(cytoline,'skel');      %removes corners
cytoline=cytoline.*cytolineskel;
%cytonuc_la=imerode(cri_la,strel('disk',2,0)); %make sure calc inside true nuc
cytoring=cro_la-cri_la;

%{
cytoringoutermass=imopen(DAs_pad,strel('disk',nucr/2,0));  %remove small debris
cytoringoutermass=bwmorph(cytoringoutermass,'thicken',outerrad);
cytoringoutermass=cytoringoutermass & zeroborder;
%}
%{
cro_la=bwlabel(cytoringoutermass);
numcells=max(max(cro_la));
crm_la=imerode(cro_la,strel('disk',1,0));
cytoline=cro_la-crm_la;
cri_la=imerode(cro_la,strel('disk',ringwidth,0));
cytonuc_la=imerode(cro_la,strel('disk',ringmargin,0)); %make sure calc inside true nuc
realnuc_la=cro_la.*DAs_pad;
cytoring=cro_la-cri_la;
%}

%{
cytoringmidmass=imerode(cytoringoutermass,strel('disk',1,0));
cytoringmidmass=~bwmorph(~cytoringmidmass,'diag');      %break connections
cytoringmidmass=~bwmorph(~cytoringmidmass,'bridge');    %break connections
cytoringmidmass=~bwmorph(cytoringmidmass,'spur');       %remove spurs
crm_la=bwlabel(cytoringmidmass);
numcells=max(max(crm_la));
cro_la=imdilate(crm_la,strel('disk',1,0));
cri_la=imerode(crm_la,strel('disk',1,0));
cytoline=crm_la-cri_la;
cytolineskel=bwmorph(cytoline,'skel');      %removes corners
cytoline=cytoline.*cytolineskel;
%}

%{
DAs_pad=DAs_pad & zeroborder;   %necessary to imerode from edges
realnuc_la=bwlabel(DAs_pad);
numcells=max(max(realnuc_la));
smoothmask=imopen(DAs_pad,strel('disk',nucr/4,0));
realnuc_la=realnuc_la.*smoothmask;  %maintains IDs, but smoothens & removes debris
cro_la=imdilate(realnuc_la,strel('disk',1,0));
cri_la=imerode(realnuc_la,strel('disk',1,0));
cytoline=realnuc_la-cri_la;
cytolineskel=bwmorph(cytoline,'skel');      %removes corners
cytoline=cytoline.*cytolineskel;
%{
if ringwidth==3
    cri_la=imerode(cri_la,strel('disk',1,0));
end
%}
cytonuc_la=imerode(cri_la,strel('disk',2,0)); %make sure calc inside true nuc
cytoring=cro_la-cri_la;

DAs_da=regionprops(cytonuc_la,'Area','Centroid','PixelIdxList'); %finds the centroid,etc of each labeled object  %to test, type "DAs_da(1).Area"
%}

%%% assign ring segment indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seglength=5;
maxsegs=1000*(nucr/8)/seglength;  %assumes 1000 is longest perimeter possible for nucr=8
rsmap=zeros(numcells,maxsegs);
segidxcur=1;    %index 1 will be NaN to ignore for cell stats
segmentedringmap=zeros(height,width);
for i=1:numcells
    %fprintf('i = %0.0f\n',i);
    if DAs_da(i).Area==0
        continue
    end
    clidx=find(cytoline==i);
    if isempty(clidx)
        continue
    end
    [r,c]=find(cytoline==i);
    set=[c,r];  %adjust to x-y convention
    [order,status]=orderperimeter(set);
    if status==0
        continue
    end
    clidx=clidx(order);
    segidx=ceil([1:length(clidx)]/seglength)+segidxcur;
    segmentedringmap(clidx)=segidx;     %set pixels to ringsegment id
    usi=unique(segidx);
    rsmap(i,:)=[usi zeros(1,maxsegs-length(usi))];
    segidxcur=segidx(end);
end
maxsegnum=find(sum(rsmap,1)==0,1)-1;
rsmap=rsmap(:,1:maxsegnum);
rsmap(rsmap==0)=1;   %these will get a mean values of NaN and be ignored in calculations

%%% clarify boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srm=imdilate(segmentedringmap,strel('disk',1,8));
srm=srm.*logical(cytoring);
%borders=bwmorph(cytoringoutermass,'bothat');
%borders=imdilate(borders,strel('disk',ringwidth,8));
srm=srm.*~borders;
srm=srm.*ringzeroborder;

REs_mask=REs_bs>0;
REs_mask=imopen(REs_mask,strel('disk',nucr,8));
srm=srm.*REs_mask;
%}
%%% select ring segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REs_bs(REs_bs<1)=1; REs_bs_log=log(REs_bs);
srm_mean=cell2mat(struct2cell(regionprops(srm,REs_bs_log,'MeanIntensity'))');
segidxrem=segidxcur-length(srm_mean);
srm_mean=[srm_mean;ones(segidxrem,1)*NaN];
ringsegmean=srm_mean(rsmap);

%%% determine outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxsegnum=size(rsmap,2);
[ringsegdev,rslb,~]=outliers(ringsegmean,maxsegnum);
rslowliers=ringsegdev<rslb;
%%% recalc median and filter outliers %%%%%%%
ringsegmean(rslowliers)=NaN;

%{
[ringsegdev,rslb,rsub]=outliers(ringsegmean,maxsegnum);
rslowliers=ringsegdev<rslb;
rsinliers=ringsegdev<=rsub;
%ringsegmean(rslowliers)=NaN;
%[ringsegdev,~,rsub]=outliers(ringsegmean,maxsegnum);
%rsinliers=ringsegdev<=rsub;
%}

[ringsegdev,rslb,rsub]=outliers(ringsegmean,maxsegnum);
rsoutliers=ringsegdev<rslb | ringsegdev>rsub;
ringsegmean(rsoutliers)=NaN;
[ringsegdev,rslb,rsub]=outliers(ringsegmean,maxsegnum);
rsinliers=ringsegdev>=rslb & ringsegdev<=rsub;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
rshighliers=ringsegdev>rsub;
rslowliers=ringsegdev<rslb;
sumh=sum(rshighliers,2); summ=sum(rsinliers,2); suml=sum(rslowliers,2);
%%%%%%% criteria %%%%%%%%%%%%%%%%%%%%%%
midgood=sumh<=1 | summ>=4;
midhighgood=sumh>=2 | summ<=2;
highgood=(sumh>=2 & suml==0) | sumh>=4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rschoice(midgood,:)=rsinliers(midgood,:);
rschoice(midhighgood,:)=rsinliers(midhighgood,:) + rshighliers(midhighgood,:);
rschoice(highgood,:)=rshighliers(highgood,:);
rschoice=rsinliers;
%}
srm_choice=rsmap(logical(rsinliers));
goodsegmask=ismember(srm,srm_choice);
finalcytoring=cytoring.*goodsegmask;

%%%%%%%% reconstitute absent rings %%%%%%%%
cru=unique(cytoring);
fcru=unique(finalcytoring);
noring=cru(~ismember(cru,fcru));
cytoring_rzb=cytoring.*ringzeroborder;
if ~isempty(noring)
    for i=noring'
        finalcytoring(cytoring_rzb==i)=i;
    end
end

%{
%%% visualization for debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tempframe=imadjust(mat2gray(bwmorph(DAs_pad,'remove')));
tempframe=imadjust(mat2gray(REs_bs));
tempframe(:,:,2)=imadjust(mat2gray(logical(finalcytoring)));
%tempframe(:,:,3)=imadjust(mat2gray(legitedges));
tempframe(:,:,3)=imadjust(mat2gray(bwmorph(DAs_pad,'remove')));
imshow(tempframe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
end

function [ringsegdev,rslb,rsub]=outliers(ringsegmean,maxsegnum)
ringsegmed=nanmedian(ringsegmean,2)*ones(1,maxsegnum);
ringsegdev=ringsegmean-ringsegmed;
ringsegstats=ringsegdev(:);
ringsegstats(isnan(ringsegstats))=[];
rsqrt=prctile(ringsegstats,[25 75]); rsiqr=iqr(ringsegstats);
rslb=rsqrt(1)-1.5*rsiqr; rsub=rsqrt(2)+1.5*rsiqr;
end