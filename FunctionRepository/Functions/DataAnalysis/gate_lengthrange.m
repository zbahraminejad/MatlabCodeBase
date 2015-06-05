function badtraces=gate_pipdeg(pipdeg,tracestats,minlength,maxthreshold,minthreshold)
%hist(min(traces2gating,[],2),0:20:2020); xlim([0 2000]);
%hist(min(traces2gating,[],2),0:80:8080); xlim([0 8000]);
%hist(log2(min(traces2gating,[],2)),0:0.1:15.1); xlim([0 15]);
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<=minlength;
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
nosignal=max(pipdeg,[],2)<=maxthreshold;
%%% remove traces where min nuc intensity is too high %%%%%%%%%%%%%%%%%%%%%
strange=min(pipdeg,[],2)>=minthreshold;

badtraces=shorttraces | nosignal | strange;