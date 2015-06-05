function [DAs_da,finalcytoring,realnuc_la]=buildcytoring_outlier(DAs_pad,REs_bs,nucr)

if nucr==16
    outerrad=4;
    ringwidth=3;
elseif nucr==8
    outerrad=3;
    ringwidth=2;
end
ringmargin=outerrad+1;

[height,width]=size(DAs_pad);
cytoringoutermass=bwmorph(DAs_pad,'thicken',outerrad);
zeroborder=ones(height,width);
zeroborder(1,:)=0;zeroborder(end,:)=0;zeroborder(:,1)=0;zeroborder(:,end)=0;
ringzeroborder=ones(height,width);
ringzeroborder(1:ringwidth+1,:)=0;ringzeroborder(end-ringwidth:end,:)=0;ringzeroborder(:,1:ringwidth+1)=0;ringzeroborder(:,end-ringwidth:end)=0;
cytoringoutermass=cytoringoutermass & zeroborder;
cro_la=bwlabel(cytoringoutermass);
cri_la=imerode(cro_la,strel('disk',ringwidth,0));     %change this back before running! or not!
cytonuc_la=imerode(cro_la,strel('disk',ringmargin,0)); %make sure calc inside true nuc
realnuc_la=cro_la.*DAs_pad;
cytoring=cro_la-cri_la;

DAs_da=regionprops(cytonuc_la,'Area','Centroid','PixelIdxList'); %finds the centroid,etc of each labeled object  %to test, type "DAs_da(1).Area"
ringxypos=regionprops(cytoring, 'PixelIdxList','PixelList');   %PixelList=[coords from left, coords from top]

%%% find spoke bases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(DAs_da,1);
tmap=zeros(height,width);
segmentedringmap=tmap;
All_home=tmap;
for i=1:numcells
    if DAs_da(i).Area==0
        continue
    end
    cenx=round(DAs_da(i).Centroid(1));  %Centroid=[coords from left edge, coords from top]
    ceny=round(DAs_da(i).Centroid(2));
    ringcoords=ringxypos(i).PixelList;
    segmentedringmap(ringxypos(i).PixelIdxList(abs(abs(ringcoords(:,1)-cenx)-2*abs(ringcoords(:,2)-ceny))<=1 | abs(2*abs(ringcoords(:,1)-cenx)-abs(ringcoords(:,2)-ceny))<=1))=1;

    NScands=ringcoords(ringcoords(:,1)==cenx,2);
    All_home([min(NScands) max(NScands)],cenx)=i;
    EWcands=ringcoords(ringcoords(:,2)==ceny,1);
    All_home(ceny,[min(EWcands) max(EWcands)])=i;
    
    NESWcands=ringcoords(ringcoords(:,1)-cenx==-1*(ringcoords(:,2)-ceny),:);
    [~,minidx]=min(NESWcands,[],1); %returns minimum of each column
    tempidx=NESWcands(minidx(2),:);
    All_home(tempidx(2),tempidx(1))=i;
    tempidx=NESWcands(minidx(1),:);
    All_home(tempidx(2),tempidx(1))=i;
    NWSEcands=ringcoords(ringcoords(:,1)-cenx==(ringcoords(:,2)-ceny),:);
    [~,minidx]=min(NWSEcands,[],1);[~,maxidx]=max(NWSEcands,[],1);
    tempidx=NWSEcands(minidx(2),:);
    All_home(tempidx(2),tempidx(1))=i;
    tempidx=NWSEcands(maxidx(1),:);
    All_home(tempidx(2),tempidx(1))=i;
end
%%% assign bases to rings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_home_rp=regionprops(All_home,'PixelIdxList');
segmentedringmap(logical(All_home))=0;
segmentedringmap=cytoring & ~segmentedringmap;
srm_la=bwlabel(segmentedringmap);
srmpoint_la=srm_la.*logical(All_home);
numsegs=max(max(srmpoint_la));
srmpoint_area=cell2mat(struct2cell(regionprops(srmpoint_la,'Area'))');
emptypoints=find(srmpoint_area==0);
numextrapoints=sum(srmpoint_area(srmpoint_area>1)-1);
numsegs=numsegs+numextrapoints+1;  %last index will point to NaN
extrapointsidx=find(srmpoint_area>1);
srmpoint_pix=regionprops(srmpoint_la,'PixelIdxList');
for i=emptypoints'
    srmpoint_pix(i).PixelIdxList=0;   %these don't matter, just maintain index for cell2mat
end
extrapix=zeros(numextrapoints,1);
epidx=1;
for i=extrapointsidx'
    nep=srmpoint_area(i)-2;
    extrapix(epidx:epidx+nep)=srmpoint_pix(i).PixelIdxList(2:end);
    srmpoint_pix(i).PixelIdxList(2:end)=[];
    epidx=epidx+1+nep;
end
srmpoint_pix=cell2mat(struct2cell(srmpoint_pix)');
srmpoint_pix=[srmpoint_pix;extrapix];
nanpix=length(srmpoint_pix)+1;
rsmap=zeros(numcells,8);
lo=zeros(numcells,1);
for i=1:numcells
    cellpix=All_home_rp(i).PixelIdxList;
    if isempty(cellpix)
        lo(i)=1;
    else
        cellmembers=find(ismember(srmpoint_pix,cellpix))';
        dummies=ones(1,8-numel(cellmembers))*nanpix;
        rsmap(i,:)=[cellmembers,dummies];
    end
end
rsmap(find(lo),:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(cytoringoutermass,'bothat');
borders=imdilate(borders,strel('disk',outerrad,8));
cytoringunbounded=cytoring.*~borders;

%%% select ring segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srm_la=srm_la.*ringzeroborder;
srm_la=srm_la.*logical(cytoringunbounded);
srm_mean=cell2mat(struct2cell(regionprops(srm_la,REs_bs,'MeanIntensity'))');
srmrem=numsegs-length(srm_mean);
if srmrem>0
    srm_mean=[srm_mean;ones(srmrem,1)*NaN];
end
ringsegmean=srm_mean(rsmap);

%%%%%%%% determine outliers %%%%%%%%%%%%%%%
ringsegmed=nanmedian(ringsegmean,2)*ones(1,8);
ringsegdev=ringsegmean-ringsegmed;
ringsegstats=ringsegdev(:);
ringsegstats(isnan(ringsegstats))=[];
rsqrt=prctile(ringsegstats,[25 75]); rsiqr=iqr(ringsegstats);
rslb=rsqrt(1)-1.5*rsiqr; rsub=rsqrt(2)+1.5*rsiqr;
rsinliers=ringsegdev>=rslb & ringsegdev<=rsub;
rshighliers=ringsegdev>rsub;
rslowliers=ringsegdev<rslb;
sumh=sum(rshighliers,2); summ=sum(rsinliers,2); suml=sum(rslowliers,2);
midgood=sumh<=1 | summ>=4;
midhighgood=sumh>=2 | summ<=2;
highgood=(sumh>=2 & suml==0) | sumh>=4;
rschoice(midgood,:)=rsinliers(midgood,:);
rschoice(midhighgood,:)=rsinliers(midhighgood,:) + rshighliers(midhighgood,:);
rschoice(highgood,:)=rshighliers(highgood,:);

srm_choice=rsmap(logical(rschoice));
goodsegmask=ismember(srm_la,srm_choice);
finalcytoring=cytoringunbounded.*goodsegmask;
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