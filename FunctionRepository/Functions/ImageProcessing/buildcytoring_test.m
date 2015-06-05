function [DAs_da,realnuc_la,finalcytoring]=buildcytoring(DAs_pad,REs_bs,nucr)
%%% define cytosolic ring radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
midrad=3;  %two pixels past the originally identified nuclear border (cuz eroded by one pixel)
ringwidth=midrad;
[height,width]=size(DAs_pad);

%%% define border margins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeroborder=ones(height,width);
zeroborder(1,:)=0;zeroborder(end,:)=0;zeroborder(:,1)=0;zeroborder(:,end)=0;
ringzeroborder=ones(height,width);
ringzeroborder(1:ringwidth+1,:)=0;ringzeroborder(end-ringwidth:end,:)=0;ringzeroborder(:,1:ringwidth+1)=0;ringzeroborder(:,end-ringwidth:end)=0;

%%% define nuclear objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad & zeroborder;   %necessary to imerode from edges

%realnuc_la=bwlabel(DAs_pad);
%smoothmask=imopen(DAs_pad,strel('disk',round(nucr/3),0));
%realnuc_la=realnuc_la.*smoothmask;  %maintains IDs, but smoothens & removes debris

DAs_pad=imopen(DAs_pad,strel('disk',round(nucr/5),0));
DAs_pad=imerode(DAs_pad,strel('disk',1,0));     %erode by one pixel for a better nuclear intensity calc
[realnuc_la,numcells]=bwlabel(DAs_pad);
DAs_da=regionprops(realnuc_la,'Area','Centroid','PixelIdxList');


%%% define midband of cytoring with same labels as nuclei %%%%%%%%%%%%%%%%%
crm_la_org=bwlabel(bwmorph(realnuc_la,'thicken',midrad));
crm_la=zeros(height,width);
for i=1:numcells
    %if DAs_da(i).Area==0    %was previously debris
    %    continue
    %end
    x=round(DAs_da(i).Centroid(1));
    y=round(DAs_da(i).Centroid(2));
    crm_label=crm_la_org(y,x);
    if crm_label>0
        crm_la(crm_la_org==crm_label)=i;
    else                    %this only would occur if the cell was split in the middle
        DAs_da(i).Area=0;   %this will effectively remove the cell from analysis
    end
end
smoothmask=imopen(crm_la,strel('disk',round(nucr/5),0));
crm_la=crm_la.*logical(smoothmask);  %maintains IDs, but smoothens

%%% define entire cytoring %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cro_la=imdilate(crm_la,strel('disk',2,0));  %outerrad two pixels greater than midrad
cri_la=imerode(crm_la,strel('disk',1,0));   %innerrad one pixel less than midrad
cytoline=crm_la-cri_la;                     %defines exactly the cytoring midline (with label ID)
cytolineskel=bwmorph(cytoline,'skel');      %removes corners
cytoline=cytoline.*cytolineskel;            %retain label ID
%crz_la=imerode(crm_la,strel('disk',2,0));   %zerorad two pixels less than midrad
%cytoring=cro_la-crz_la;                     %define cytoring with label ID
cytoring=cro_la-cri_la;                     %define cytoring with label ID.  3 pixels wide

%%% assign ring segment indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seglength=round(nucr/2);
maxsegs=round(1000*nucr/seglength);  %give enough room for many segments
rsmap=zeros(numcells,maxsegs);
segidxcur=1;                         %index 1 will be NaN to ignore for cell stats
segmentedringmap=zeros(height,width);
for i=1:numcells
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

%%% detect cell-cell borders and inflate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(crm_la_org,'bothat');
borders=imdilate(borders,strel('disk',ringwidth,8));

%%% clarify boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srm=imdilate(segmentedringmap,strel('disk',1,8));
srm=srm.*logical(cytoring);
srm=srm.*~borders;
srm=srm.*ringzeroborder;

%%% select ring segments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REs_bs(REs_bs<1)=1; REs_bs_log=log(REs_bs);
srm_mean=cell2mat(struct2cell(regionprops(srm,REs_bs,'MeanIntensity'))');
segidxrem=segidxcur-length(srm_mean);
srm_mean=[srm_mean;ones(segidxrem,1)*NaN];  %correct indexing retained, of course
ringsegmean=srm_mean(rsmap);

%%% determine outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxsegnum=size(rsmap,2);
[ringsegdev,rslb,~]=outliers(ringsegmean,maxsegnum);
rsinliers=ringsegdev>=rslb;
srm_choice=rsmap(logical(rsinliers));                   %remove ouliers
goodsegmask=ismember(srm,srm_choice);
finalcytoring=cytoring.*goodsegmask;

%%%%%%%% reconstitute absent rings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%% visualization for debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tempframe=imadjust(mat2gray(bwmorph(DAs_pad,'remove')));
tempframe=imadjust(mat2gray(REs_bs));
tempframe(:,:,2)=imadjust(mat2gray(logical(finalcytoring)));
%tempframe(:,:,3)=imadjust(mat2gray(legitedges));
tempframe(:,:,3)=imadjust(mat2gray(bwmorph(DAs_pad,'remove')));
imshow(tempframe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
end

function [ringsegdev,rslb,rsub]=outliers(ringsegmean,maxsegnum)
ringsegmed=nanmedian(ringsegmean,2)*ones(1,maxsegnum);
ringsegdev=ringsegmean-ringsegmed;
ringsegstats=ringsegdev(:);
ringsegstats(isnan(ringsegstats))=[];
rsqrt=prctile(ringsegstats,[25 75]); rsiqr=iqr(ringsegstats);
rslb=rsqrt(1)-0.5*rsiqr; rsub=rsqrt(2)+1.5*rsiqr;
end