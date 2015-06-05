function mask=segmentdeflections_2(mask,nucr,cellsizeratio,debrisarea)
bigthresh=round(cellsizeratio*pi*nucr^2);  %1.33*pi*nucr^2 is the avg cell size
bigmask=bwareaopen(mask,bigthresh);

[B,L]=bwboundaries(bigmask,'noholes');
obnum=max(L(:));
bordermask=zeros(size(mask));
for ci=1:obnum
    scmask=L==ci;
    bordermask=splitdeflections_2(B{ci},bordermask,nucr,scmask);
end

mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');
mask=bwareaopen(mask,debrisarea); %bin1:0.75 bin2:0.5
end