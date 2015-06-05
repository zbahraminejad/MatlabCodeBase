function [DAs_da,finalcytoring,realnuc_la]=buildcytoring(DAs_pad,REs_bs,nucr)

if nucr==16
    outerrad=4;
    ringwidth=3;
elseif nucr==8
    outerrad=3;
    ringwidth=2;
end
ringmargin=outerrad+1;
[height,width]=size(DAs_pad);
%REs_ma=ThreshImage_MC(REs_bs,0);  %looser threshold
REs_ma=ThreshImage_MC(REs_bs,5);  %optimal threshold
REs_ma=imopen(REs_ma,strel('disk',nucr/4,8));
REs_ma=bwmorph(REs_ma,'remove');
REs_ma(REs_ma & DAs_pad)=0;
%outerrad=nucr/4;
%ringwidth=outerrad/2+1;
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

%%% find spoke bases %%%%%%%%
numcells=size(DAs_da,1);
tmap=zeros(height,width);
segmentedringmap=tmap;
N_home=tmap;NE_home=tmap;E_home=tmap;SE_home=tmap;S_home=tmap;SW_home=tmap;W_home=tmap;NW_home=tmap;
for i=1:numcells
    if DAs_da(i).Area==0
        continue
    end
    cenx=round(DAs_da(i).Centroid(1));  %Centroid=[coords from left edge, coords from top]
    ceny=round(DAs_da(i).Centroid(2));
    ringcoords=ringxypos(i).PixelList;
    segmentedringmap(ringxypos(i).PixelIdxList(abs(abs(ringcoords(:,1)-cenx)-2*abs(ringcoords(:,2)-ceny))<=1 | abs(2*abs(ringcoords(:,1)-cenx)-abs(ringcoords(:,2)-ceny))<=1))=1;

    NScands=ringcoords(ringcoords(:,1)==cenx,2); 
    N_home(min(NScands),cenx)=i; S_home(max(NScands),cenx)=i;
    EWcands=ringcoords(ringcoords(:,2)==ceny,1); 
    E_home(ceny,max(EWcands))=i; W_home(ceny,min(EWcands))=i;

    NESWcands=ringcoords(ringcoords(:,1)-cenx==-1*(ringcoords(:,2)-ceny),:);
    [~,minidx]=min(NESWcands,[],1); %returns minimum of each column
    tempidx=NESWcands(minidx(2),:);
    NE_home(tempidx(2),tempidx(1))=i; 
    tempidx=NESWcands(minidx(1),:);
    SW_home(tempidx(2),tempidx(1))=i;
    NWSEcands=ringcoords(ringcoords(:,1)-cenx==(ringcoords(:,2)-ceny),:);
    [~,minidx]=min(NWSEcands,[],1);[~,maxidx]=max(NWSEcands,[],1);
    tempidx=NWSEcands(minidx(2),:);
    NW_home(tempidx(2),tempidx(1))=i;
    tempidx=NWSEcands(maxidx(1),:);
    SE_home(tempidx(2),tempidx(1))=i;
end
All_home=N_home + NE_home + E_home + SE_home + S_home + SW_home + W_home + NW_home;
%segmentedringmap=imdilate(segmentedringmap,strel('disk',1,8));
segmentedringmap=cytoring & ~segmentedringmap;
%segmentedringmap=imopen(segmentedringmap,strel('disk',1,8));  %was left in for first 20x
srm_la=bwlabel(segmentedringmap);
%%% establish cytosol gradients %%%%%%%%%%%%%%%%%%%
%{
cytogradients=edge(REs_bs,'canny');
cytogradients=bwmorph(cytogradients,'thin',Inf);
cytogradients=bwmorph(cytogradients,'diag'); %removes 8-connectivity of background
cg_da=regionprops(cytogradients,'BoundingBox','PixelIdxList');
for i=1:size(cg_da,1)
    longestbound=max(cg_da(i).BoundingBox([3 4]));
    if longestbound < nucr*2
        cytogradients(cg_da(i).PixelIdxList)=0;
    end
end
%}
%%% remove ring intersected by translocating sensor gradients %%%%%%%%%
boundthresh=nucr/4;
%{
legitedges=bwlabel(cytogradients);
gradientratio=3;
innerscreen=legitedges.*cytoringoutermass;
outerscreen=legitedges.*~cytoringoutermass;
pseudonuc=regionprops(innerscreen,'Area');
pseudocyto=regionprops(outerscreen,'Area');
highestinner=max(max(innerscreen));
highestouter=max(max(outerscreen));
for j=1:highestinner
    if j>highestouter
        legitedges(legitedges==j)=0;
    elseif pseudonuc(j).Area > gradientratio*pseudocyto(j).Area
        legitedges(legitedges==j)=0;
    end
end

legitedges=legitedges | REs_ma;
inbounds=cytoringoutermass & ~legitedges;
%}
inbounds=cytoringoutermass & ~REs_ma;
inbounds=imopen(inbounds,strel('disk',boundthresh,8));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(cytoringoutermass,'bothat');
borders=imdilate(borders,strel('disk',outerrad,8));
cytoringunbounded=cytoring.*~borders;
cytoringbounded=cytoringunbounded.*inbounds;

%%% determine valid ring segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_potential=cytoringunbounded & All_home;
potentialsegments=unique(srm_la.*All_potential);
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

N_exm=detectoverlaps(N_home,N_out,1,cytoringoutermass,REs_ma,REs_bs);
NW_exm=detectoverlaps(NW_home,NW_out,1,cytoringoutermass,REs_ma,REs_bs);
W_exm=detectoverlaps(W_home,W_out,1,cytoringoutermass,REs_ma,REs_bs);
SW_exm=detectoverlaps(SW_home,SW_out,1,cytoringoutermass,REs_ma,REs_bs);
S_exm=detectoverlaps(S_home,S_out,-1,cytoringoutermass,REs_ma,REs_bs);
SE_exm=detectoverlaps(SE_home,SE_out,-1,cytoringoutermass,REs_ma,REs_bs);
E_exm=detectoverlaps(E_home,E_out,-1,cytoringoutermass,REs_ma,REs_bs);
NE_exm=detectoverlaps(NE_home,NE_out,-1,cytoringoutermass,REs_ma,REs_bs);

excludemap=N_exm | NW_exm | W_exm | SW_exm | S_exm | SE_exm | E_exm | NE_exm;

excludesegments=unique(srm_la.*excludemap);
potentialsegments(ismember(potentialsegments,excludesegments))=[];
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


%%% visualization for debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tempframe=imadjust(mat2gray(logical(realnuc_la)));
tempframe=imadjust(mat2gray(REs_bs));
tempframe(:,:,2)=imadjust(mat2gray(logical(finalcytoring)));
tempframe(:,:,3)=imadjust(mat2gray(REs_ma));
imshow(tempframe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
end