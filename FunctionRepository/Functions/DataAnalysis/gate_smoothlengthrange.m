function [signals,badtraces]=gate_smoothlengthrange(tracedata,tracestats,channel,minlength,maxthreshold,minthreshold)
%hist(min(signals,[],2),0:20:2020); xlim([0 2000]);
%hist(max(signals,[],2),0:10:1010); xlim([0 1000]);
%hist(log2(min(signals,[],2)),0:0.1:15.1); xlim([0 15]);

%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
signals=tracedata(:,:,channel);
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<=minlength;
nosignal=max(signals,[],2)<=maxthreshold;
strange=min(signals,[],2)>=minthreshold;
badtraces=shorttraces | nosignal | strange;