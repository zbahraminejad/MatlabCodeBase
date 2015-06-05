function mask=getnucmask_status(DAs_bs,nucr)
%%% Perform nuclear segmentation by 'Concave Detection' method %%%%%%%%%%%%
% Author: Mingyu Chung
% Last Revision: 2/5/2013
%
% How to use this script:
% Arg1: background-subtracted image of cell nuclei (e.g. Hoechst stain)
% Arg2: average nuclear radius in pixels. This determines the minimum size
%       of nuclei to not consider as debris, and it only slightly matters
%       for segmenting the nuclei.
%
% What this script does:
% 1. Generates a mask purely based on intensity threshold
% 2. Cleans up the mask by breaking some strands, spurs, etc.
% 3. For each object, searches perimeter and detects concave deflections
% 4. Bridges deflections/vertices to form segmentation borders
%    *Criteria applied for multiple vertices (e.g. more than 2 cells
%     overlapping, etc.)
% 
% This script requires:
% - ThreshImage_MC.m
% - orderperimeter.m
% - segmentnuclei.m
% - bridge.m
%
% Users:
% Gautam Dey (gdey@stanford.edu)
% Anshul Rana (anshulr@stanford.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% define initial mask by threshold intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=ThreshImage_MC(DAs_bs,0);            %2nd arg = # bins above threshold
mask=imfill(mask,'holes');
mask=padarray(mask,[1 1]);                %necessary for imopen
mask=imopen(mask,strel('disk',nucr/2,0)); %this would remove debris & spurs
mask=~bwmorph(~mask,'diag');              %break diagonal connections
mask=~bwmorph(~mask,'bridge');            %break orthogonal connections
nucedges=bwmorph(mask,'remove');          %remove interior pixels
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(mask));

%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];                  %adjust to x-y convention
    [order,status]=orderperimeter(coorset);
    if status==0                    %error, skip segmentation for this cell
        continue
    end
    orderedset=coorset(order,:);
    bordermask=segmentnuclei(orderedset,bordermask,nucr);
end

%%% Incorporate segment borders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=mask(2:end-1,2:end-1);
bordermask=bordermask(2:end-1,2:end-1);
mask=mask & ~bordermask;                  %add segmentations
mask=~bwmorph(~mask,'diag');              %remove diagonal connections
end