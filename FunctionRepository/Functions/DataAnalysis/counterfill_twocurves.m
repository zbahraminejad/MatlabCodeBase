function counterfill(bin,topvals,bottomvals,bstep,lowbin,highbin,colorchar)
hold on;
[~,idxlow]=min(abs(bin-lowbin));
[~,idxhigh]=min(abs(bin-highbin));
midlow=bin(idxlow); midhigh=bin(idxhigh);
midbin=midlow:bstep:midhigh;
midfill=[midbin fliplr(midbin)];
fill(midfill,[topvals(idxlow:idxhigh);bottomvals(idxhigh:-1:idxlow)],colorchar,'FaceAlpha',0.7);