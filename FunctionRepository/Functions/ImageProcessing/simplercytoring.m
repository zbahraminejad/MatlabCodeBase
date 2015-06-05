function [DAs_da,realnuc_la,finalcytoring]=simplercytoring(DAs_pad,nucr)
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

%%% define midband of cytoring with same labels as nuclei %%%%%%%%%%%%%%%%%
ringedge=imdilate(realnuc_la,strel('disk',3,0));  %outerrad one pixel greater than midrad
cytoring=ringedge-realnuc_la;                     %define cytoring with label ID

%%% detect cell-cell borders and inflate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
borders=bwmorph(DAs_pad,'bothat');
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
tempframe(:,:,3)=imadjust(mat2gray(bwmorph(DAs_pad,'remove')));
imshow(tempframe);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
end