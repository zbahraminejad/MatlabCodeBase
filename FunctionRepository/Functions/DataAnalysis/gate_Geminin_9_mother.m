function [signals,badtraces]=gate_Geminin_9_mother(tracedata,tracestats,motherstats,channel,maxthreshold,minthreshold,maxpos,minpos)
%hist(max(signals,[],2),100);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
signals=tracedata(:,:,channel);
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toolow=ones(numtraces,1)*NaN;
switch maxpos
    case 0
        toolow=max(signals,[],2)<maxthreshold;
    case 1
        for i=1:numtraces
            %toolow(i)=signals(i,tracestats(i,1))<maxthreshold;
            toolow(i)=max(signals(i,motherstats(i,1):motherstats(i,2)))<maxthreshold;
        end
    case 2
        for i=1:numtraces
            toolow(i)=signals(i,tracestats(i,2))<maxthreshold;
        end
end

toohigh=ones(numtraces,1)*NaN;
switch minpos
    case 0
        toohigh=min(signals,[],2)>minthreshold;
    case 1
        for i=1:numtraces
            toohigh(i)=min(signals(i,motherstats(i,1):motherstats(i,2)))>minthreshold;
        end
    case 2
        for i=1:numtraces
            toohigh(i)=min(signals(i,tracestats(i,1):tracestats(i,2)))>minthreshold;
        end
end

badtraces=toolow | toohigh;
signals=tracedata(:,:,channel); %return raw signal (not smoothened)