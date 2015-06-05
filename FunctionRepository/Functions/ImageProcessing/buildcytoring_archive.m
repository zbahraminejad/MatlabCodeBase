function [DAs_da,finalcytoring,realnuc_la]=buildcytoring_outlier_test(DAs_pad,REs_bs,nucr)

if nucr==16
    outerrad=4;
    ringwidth=3;
elseif nucr==8
    outerrad=3;
    ringwidth=2;
end
ringmargin=outerrad+1;

[height,width]=size(DAs_pad);
%{
%REs_ma=ThreshImage_MC(REs_bs,0);  %looser threshold
REs_ma=ThreshImage_MC(REs_bs,5);  %optimal threshold
%REs_ma=imopen(REs_ma,strel('disk',nucr/4,8));
REs_ma=imopen(REs_ma,strel('disk',nucr,8));   %normally nucr/4
REs_ma=bwmorph(REs_ma,'remove');
REs_ma(REs_ma & DAs_pad)=0;
%}
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
%DAs_da=regionprops(realnuc_la,'Area','Centroid','PixelIdxList');
ringxypos=regionprops(cytoring, 'PixelIdxList','PixelList');   %PixelList=[coords from left, coords from top]

%%% find spoke bases %%%%%%%%
numcells=size(DAs_da,1);
tmap=zeros(height,width);
segmentedringmap=tmap;
All_home=tmap;
%N_home=tmap;NE_home=tmap;E_home=tmap;SE_home=tmap;S_home=tmap;SW_home=tmap;W_home=tmap;NW_home=tmap;
for i=1:numcells
    if DAs_da(i).Area==0
        continue
    end
    cenx=round(DAs_da(i).Centroid(1));  %Centroid=[coords from left edge, coords from top]
    ceny=round(DAs_da(i).Centroid(2));
    ringcoords=ringxypos(i).PixelList;
    segmentedringmap(ringxypos(i).PixelIdxList(abs(abs(ringcoords(:,1)-cenx)-2*abs(ringcoords(:,2)-ceny))<=1 | abs(2*abs(ringcoords(:,1)-cenx)-abs(ringcoords(:,2)-ceny))<=1))=1;

    NScands=ringcoords(ringcoords(:,1)==cenx,2);
    %N_home(min(NScands),cenx)=i; S_home(max(NScands),cenx)=i;
    %temp_home(1,:)=[min(NScands),cenx];
    %temp_home(5,:)=[max(NScands),cenx);
    All_home([min(NScands) max(NScands)],cenx)=i;
    EWcands=ringcoords(ringcoords(:,2)==ceny,1); 
    %E_home(ceny,max(EWcands))=i; W_home(ceny,min(EWcands))=i;
    %temp_home(3,:)=[ceny,max(EWcands)];
    %temp_home(7,:)=[ceny,min(EWcands)];
    All_home(ceny,[min(EWcands) max(EWcands)])=i;
    
    NESWcands=ringcoords(ringcoords(:,1)-cenx==-1*(ringcoords(:,2)-ceny),:);
    [~,minidx]=min(NESWcands,[],1); %returns minimum of each column
    tempidx=NESWcands(minidx(2),:);
    %NE_home(tempidx(2),tempidx(1))=i;
    %temp_home(2,:)=[tempidx(2),tempidx(1)];
    All_home(tempidx(2),tempidx(1))=i;
    tempidx=NESWcands(minidx(1),:);
    %SW_home(tempidx(2),tempidx(1))=i;
    %temp_home(6,:)=[tempidx(2),tempidx(1)];
    All_home(tempidx(2),tempidx(1))=i;
    NWSEcands=ringcoords(ringcoords(:,1)-cenx==(ringcoords(:,2)-ceny),:);
    [~,minidx]=min(NWSEcands,[],1);[~,maxidx]=max(NWSEcands,[],1);
    tempidx=NWSEcands(minidx(2),:);
    %NW_home(tempidx(2),tempidx(1))=i;
    %temp_home(8,:)=[tempidx(2),tempidx(1)];
    All_home(tempidx(2),tempidx(1))=i;
    tempidx=NWSEcands(maxidx(1),:);
    %SE_home(tempidx(2),tempidx(1))=i;
    %temp_home(4,:)=[tempidx(2),tempidx(1)];
    All_home(tempidx(2),tempidx(1))=i;
end
%All_home=N_home + NE_home + E_home + SE_home + S_home + SW_home + W_home + NW_home;
%All_home_rp=cell2mat(struct2cell(regionprops(All_home,'PixelIdxList')));
All_home_rp=regionprops(All_home,'PixelIdxList');
%segmentedringmap=imdilate(segmentedringmap,strel('disk',1,8));
segmentedringmap(logical(All_home))=0;
segmentedringmap=cytoring & ~segmentedringmap;
%segmentedringmap=imopen(segmentedringmap,strel('disk',1,8));
srm_la=bwlabel(segmentedringmap);
%numsegs=max(max(srm_la));
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

%{
for i=1:length(srmpoint_pix)
    if isempty(srmpoint_pix(i).PixelIdxList)\
        srmpoint_pix(i).PixelIdxList=0;
    end
end
%}
%{
sc=struct2cell(srmpoint_rp)';
scmat=zeros(length(sc),1);
for i=1:length(sc)
    fprintf('i = %0.0f\n',i);
    scmat(i)=sc{i};
end
%}
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
        %rsmap(i,:)=find(ismember(srmpoint_pix,cellpix))';
    end
end
rsmap(find(lo),:)=[];

%%% remove ring intersected by translocating sensor threshold %%%%%%%%%
%boundthresh=nucr/4;
%boundthresh=outerrad;
%inbounds=cytoringoutermass & ~REs_ma;
%inbounds=imopen(inbounds,strel('disk',boundthresh,8));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(cytoringoutermass,'bothat');
borders=imdilate(borders,strel('disk',outerrad,8));
cytoringunbounded=cytoring.*~borders;
%cytoringbounded=cytoringunbounded.*inbounds;

%{
%%% determine valid ring segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All_potential=cytoringunbounded & All_home;
%potentialsegments=unique(srm_la.*All_potential);
outspokelength=nucr*2; diagspoke=round(outspokelength*0.7);

N_home=N_home.*logical(cytoringunbounded); S_home=S_home.*logical(cytoringunbounded);
E_home=E_home.*logical(cytoringunbounded); W_home=W_home.*logical(cytoringunbounded);
NE_home=NE_home.*logical(cytoringunbounded); SW_home=SW_home.*logical(cytoringunbounded);
NW_home=NW_home.*logical(cytoringunbounded); SE_home=SE_home.*logical(cytoringunbounded);

base=zeros(2*outspokelength+1,1);
nhood=base; nhood(1:outspokelength)=1; N_out=imdilate(N_home,strel(nhood)); N_out(N_out<0)=0;
nhood=base; nhood(outspokelength+2:end)=1; S_out=imdilate(S_home,strel(nhood));S_out(S_out<0)=0;
base=zeros(1,2*outspokelength+1);
nhood=base; nhood(outspokelength+2:end)=1; E_out=imdilate(E_home,strel(nhood));E_out(E_out<0)=0;
nhood=base; nhood(1:outspokelength)=1; W_out=imdilate(W_home,strel(nhood));W_out(W_out<0)=0;
base=eye(2*diagspoke+1,2*diagspoke+1);
nhood=base; nhood(diagspoke+1:end,diagspoke+1:end)=0; NW_out=imdilate(NW_home,strel(nhood));NW_out(NW_out<0)=0;
nhood=base; nhood(1:diagspoke+1,1:diagspoke+1)=0; SE_out=imdilate(SE_home,strel(nhood));SE_out(SE_out<0)=0;
base=flipud(base);
nhood=base; nhood(diagspoke+1:end,1:diagspoke+1)=0; NE_out=imdilate(NE_home,strel(nhood));NE_out(NE_out<0)=0;
nhood=base; nhood(1:diagspoke+1,diagspoke+1:end)=0; SW_out=imdilate(SW_home,strel(nhood));SW_out(SW_out<0)=0;

All_out=N_out | NE_out | E_out | SE_out | S_out | SW_out | W_out | NW_out;

N_exm=detectoverlaps(N_home,N_out,1,cytoringoutermass,REs_ma,REs_bs,height);
NW_exm=detectoverlaps(NW_home,NW_out,1,cytoringoutermass,REs_ma,REs_bs,height);
W_exm=detectoverlaps(W_home,W_out,1,cytoringoutermass,REs_ma,REs_bs,height);
SW_exm=detectoverlaps(SW_home,SW_out,1,cytoringoutermass,REs_ma,REs_bs,height);
S_exm=detectoverlaps(S_home,S_out,-1,cytoringoutermass,REs_ma,REs_bs,height);
SE_exm=detectoverlaps(SE_home,SE_out,-1,cytoringoutermass,REs_ma,REs_bs,height);
E_exm=detectoverlaps(E_home,E_out,-1,cytoringoutermass,REs_ma,REs_bs,height);
NE_exm=detectoverlaps(NE_home,NE_out,-1,cytoringoutermass,REs_ma,REs_bs,height);

overlapmap=N_exm | NW_exm | W_exm | SW_exm | S_exm | SE_exm | E_exm | NE_exm;
%}
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
%%% recalc median and filter outliers %%%%%%%
%rsoutliers=ringsegdev<rslb | ringsegdev>rsub;
%ringsegmean(rsoutliers)=NaN;
%ringsegdev=ringsegmean-ringsegmed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
rsinliers=ringsegdev>=rslb & ringsegdev<=rsub;
rsoutliers=ringsegdev<rslb | ringsegdev>rsub;
bias=sum(rsinliers,2)>=sum(rsoutliers,2);
rschoice=zeros(size(rsmap));
rschoice(bias,:)=rsinliers(bias,:);
rschoice(~bias,:)=rsoutliers(~bias,:);
%}

rsinliers=ringsegdev>=rslb & ringsegdev<=rsub;
rshighliers=ringsegdev>rsub;
rslowliers=ringsegdev<rslb;
%{
bias=sum(rshighliers,2)<3;
rschoice(bias,:)=rsinliers(bias,:);
rschoice(~bias,:)=rshighliers(~bias,:);
%}
sumh=sum(rshighliers,2); summ=sum(rsinliers,2); suml=sum(rslowliers,2);
%midgood=(sumh<=1 | (sumh==2&suml<=1)) & summ>2;
midgood=sumh<=1 | summ>=4;
%midhighgood=(sumh>=2 & suml>=2) | summ<=2;
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

%threshsegments=unique(srm_la.*REs_ma);          %find segments that cross hDHB threshold
%overlapsegments=unique(srm_la.*overlapmap);     %find segments from prior section
%blacklist=unique([threshsegments;overlapsegments]);

%{
potentialsegments(ismember(potentialsegments,overlapsegments))=[];
nonoverlapsegments=ismember(srm_la,potentialsegments);
nonoverlapsegments=nonoverlapsegments & ringzeroborder;
nonoverlap_unbounded=cytoringunbounded.*nonoverlapsegments;     %ring besides cell-cell jcn
nonoverlap_bounded=cytoringbounded.*nonoverlapsegments;                %ring besides cell-cell jcn & hDHB gradient
nonoverlap_difference=nonoverlap_unbounded.*logical(~nonoverlap_bounded);
ringsample=regionprops(nonoverlap_bounded,'Area');
minringsample=nucr*2*pi;
ringinsufficient=zeros(length(ringsample),1);
for k=1:length(ringsample)
    if ringsample(k).Area <= minringsample
        ringinsufficient(k)=1;
    end
end
ringinsufficient=find(ringinsufficient);
nonoverlap_difference=nonoverlap_difference.*ismember(nonoverlap_difference,ringinsufficient);
finalcytoring=nonoverlap_bounded + nonoverlap_difference;
%%% For any cells that have no finalcytoring, replace with original ring %%
ringexists=unique(finalcytoring);
ringexists(1)=[];   %remove zero
ringpreexists=[1:size(ringxypos,1)]';
ringmissing=ringpreexists(~ismember(ringpreexists,ringexists));
missingcytoring=cytoring.*ismember(cytoring,ringmissing);
finalcytoring=finalcytoring + missingcytoring;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

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