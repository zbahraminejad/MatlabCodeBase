function [bordermask,borderflag]=attemptsplit(i,nuc_label_mask,bordermask,nucr)
%%% attempt to segment deflections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r,c]=find(nuc_label_mask==i);
coorset=[c,r];  %adjust to x-y convention
[order,status]=orderperimeter([c,r]);
if status==0    %unable to order perimeter
    %fprintf('unable to order perimeter\n');
    borderflag=0;
    return;
end
orderedset=coorset(order,:);
[bordermask,borderflag]=splitdeflections(orderedset,bordermask,nucr);
end