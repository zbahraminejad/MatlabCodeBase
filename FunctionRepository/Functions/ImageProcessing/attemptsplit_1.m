function [bordermask,borderflag]=attemptsplit_1(i,nuc_label,bordermask,nucr)
singlecell=nuc_label==i;
B=bwboundaries(singlecell,'noholes');
[bordermask,borderflag]=splitdeflections_1(B{1},bordermask,nucr);
end