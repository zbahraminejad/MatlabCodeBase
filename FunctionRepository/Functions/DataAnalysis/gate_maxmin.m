function [signals,badtraces]=gate_maxmin(tracedata,tracestats,channel,maxpos,maxthreshold,minpos,minthreshold)
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
nosignalmin=ones(numtraces,1)*NaN;
nosignalmax=ones(numtraces,1)*NaN;
switch minpos
    case 0
        nosignalmin=min(signals,[],2)>minthreshold;
    case 1
        for i=1:numtraces
            nosignalmin(i)=signals(i,tracestats(i,1))>minthreshold;
        end
    case 2
        for i=1:numtraces
            nosignalmin(i)=signals(i,tracestats(i,2))>minthreshold;
        end
end
switch maxpos
    case 0
        nosignalmax=max(signals,[],2)<maxthreshold;
    case 1
        for i=1:numtraces
            nosignalmax(i)=signals(i,tracestats(i,1))<maxthreshold;
        end
    case 2
        for i=1:numtraces
            nosignalmax(i)=signals(i,tracestats(i,2))<maxthreshold;
        end
end

badtraces=nosignalmin | nosignalmax;
signals=tracedata(:,:,channel); %return raw signal (not smoothened)