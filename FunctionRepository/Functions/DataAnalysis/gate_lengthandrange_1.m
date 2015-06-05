function [signals,badtraces]=gate_smoothlengthrange(tracedata,tracestats,channel,minlength,maxthreshold,minthreshold,noisethresh)
%hist(min(signals,[],2),0:10:1001); xlim([0 1000]);
%hist(max(signals,[],2),0:10:2010); xlim([0 2000]);
%hist(log2(min(signals,[],2)),0:0.1:15.1); xlim([0 15]);
%hist(noise,500);

%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
signals=tracedata(:,:,channel);
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<minlength;
nosignal=max(signals,[],2)<maxthreshold;
strange=min(signals,[],2)>minthreshold;
noise=max(diff(signals,1,2),[],2);
noisy=noise>noisethresh;


badtraces=shorttraces | nosignal | strange | noisy;
signals=tracedata(:,:,channel); %return raw signal (not smoothened)