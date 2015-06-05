function badtraces=gate_pipdeg(pipdeg,tracestats,minlength,maxthreshold,minthreshold)
numtraces=size(pipdeg,1);
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<=minlength;
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
nosignal=max(pipdeg,[],2)<=maxthreshold;
%%% remove traces where min nuc intensity is too high %%%%%%%%%%%%%%%%%%%%%
strange=min(pipdeg,[],2)>=minthreshold;

badtraces=shorttraces | nosignal | strange;