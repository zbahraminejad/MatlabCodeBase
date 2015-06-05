function [signals,badtraces]=gate_pipsensor(tracedata,tracestats,channel,firstthreshold,maxthreshold,minthreshold,negthreshold)
%{
hist(max(signals,[],2),0:2:202); xlim([0 200]);
hist(max(signals,[],2),0:10:2010); xlim([0 2000]);
hist(log2(min(signals,[],2)),0:0.1:15.1); xlim([0 15]);

firstvals=ones(numtraces,1)*NaN;
for i=1:numtraces
    firstvals(i)=signals(i,tracestats(i,1));
end
hist(firstvals,100);
%}

%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
signals=tracedata(:,:,channel);
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalmin=min(signals,[],2)<minthreshold & min(signals,[],2)>negthreshold;
totalmax=max(signals,[],2)>maxthreshold;
firstmax=ones(numtraces,1)*NaN;
for i=1:numtraces
    firstmax(i)=signals(i,tracestats(i,1))>firstthreshold;
end
goodtraces=totalmin & totalmax & firstmax;
badtraces=~goodtraces;
signals=tracedata(:,:,channel); %return raw signal (not smoothened)