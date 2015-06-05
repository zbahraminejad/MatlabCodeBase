function badtraces=remove_range_noisy(signals,maxthreshold,minthreshold)
numtraces=size(signals,1);
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
nosignal=max(signals,[],2)<=maxthreshold;
%%% remove traces where min nuc intensity is too high %%%%%%%%%%%%%%%%%%%%%
strange=min(signals,[],2)>=minthreshold;
%%% remove noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numtraces
    maxval=max(signals(i,:));
    signals(i,:)=signals(i,:)/maxval;
end
swingthresh=1; %default 1
noisy=max(abs(diff(signals,2,2)),[],2)>swingthresh;

badtraces=nosignal | strange | noisy;