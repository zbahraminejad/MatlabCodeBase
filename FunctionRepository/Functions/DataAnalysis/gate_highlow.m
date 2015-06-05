function badtraces=gate_pipdeg(pipdeg,tracestats,minlength,maxthreshold,minthreshold)
numtraces=size(pipdeg,1);
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<=minlength;
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
nosignal=max(pipdeg,[],2)<=maxthreshold;
%%% remove traces where min nuc intensity is too high %%%%%%%%%%%%%%%%%%%%%
strange=min(pipdeg,[],2)>=minthreshold;
%%% remove noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numtraces
    maxval=max(pipdeg(i,:));
    pipdeg(i,:)=pipdeg(i,:)/maxval;
end
swingthresh=1; %default 1
noisy=max(abs(diff(pipdeg,2,2)),[],2)>swingthresh;

badtraces=shorttraces | nosignal | strange | noisy;