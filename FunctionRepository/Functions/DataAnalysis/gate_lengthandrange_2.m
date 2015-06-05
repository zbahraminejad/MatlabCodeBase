function [signals,badtraces]=gate_lengthandrange_2(tracedata,tracestats,channel,minlength,startrange)
%hist(min(signals,[],2),0:10:1001); xlim([0 1000]);
%hist(max(signals,[],2),0:10:2010); xlim([0 2000]);
%hist(log2(min(signals,[],2)),0:0.1:15.1); xlim([0 15]);

%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
signals=tracedata(:,:,channel);
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<minlength;
toolow=mean(signals(:,1:5),2)<startrange(1);
toohigh=mean(signals(:,1:5),2)>startrange(2);
badtraces=shorttraces | toolow | toohigh;
signals=tracedata(:,:,channel); %return raw signal (not smoothened)