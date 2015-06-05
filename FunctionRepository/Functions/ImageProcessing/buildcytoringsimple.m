function [DAs_da,realnuc_la,finalcytoring]=simplecytoring(DAs_pad,nucr)
%%% define cytosolic ring radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
midrad=2;
ringwidth=midrad+1;
[height,width]=size(DAs_pad);

%%% define border margins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeroborder=ones(height,width);
zeroborder(1,:)=0;zeroborder(end,:)=0;zeroborder(:,1)=0;zeroborder(:,end)=0;
ringzeroborder=ones(height,width);
ringzeroborder(1:ringwidth+1,:)=0;ringzeroborder(end-ringwidth:end,:)=0;ringzeroborder(:,1:ringwidth+1)=0;ringzeroborder(:,end-ringwidth:end)=0;

%%% define nuclear objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad & zeroborder;   %necessary to imerode from edges
DAs_pad=imopen(DAs_pad,strel('disk',round(nucr/5),0));
realnuc_la=bwlabel(DAs_pad);

DAs_da=regionprops(realnuc_la,'Area','Centroid','PixelIdxList');
numcells=length(DAs_da);
%%% define midband of cytoring with same labels as nuclei %%%%%%%%%%%%%%%%%
crm_la_org=bwlabel(bwmorph(realnuc_la,'thicken',midrad));
crm_la=zeros(height,width);
for i=1:numcells
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
cro_la=imdilate(crm_la,strel('disk',1,0));  %outerrad one pixel greater than midrad
crz_la=imerode(crm_la,strel('disk',2,0));   %realnuc_la could be different
cytoring=cro_la-crz_la;                     %define cytoring with label ID

%%% detect cell-cell borders and inflate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(crm_la_org,'bothat');
borders=imdilate(borders,strel('disk',ringwidth,8));

%%% clarify boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finalcytoring=cytoring.*~borders;
finalcytoring=finalcytoring.*ringzeroborder;

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