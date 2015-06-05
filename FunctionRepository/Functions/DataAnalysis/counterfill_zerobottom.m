function counterfill(bin,pdfvals,bstep,lowbin,highbin,colorchar)
hold on;
[~,idxlow]=min(abs(bin-lowbin));
[~,idxhigh]=min(abs(bin-highbin));
midlow=bin(idxlow); midhigh=bin(idxhigh);
midbin=midlow:bstep:midhigh;
midfill=[midbin fliplr(midbin)];
fill(midfill,[pdfvals(idxlow:idxhigh);zeros(numel(midbin),1)],colorchar,'FaceAlpha',0.3);