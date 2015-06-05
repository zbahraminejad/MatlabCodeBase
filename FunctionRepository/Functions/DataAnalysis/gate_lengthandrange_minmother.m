function [signals,badtraces]=gate_lengthandrange_maxpos(tracedata,tracestats,motherstats,channel,minlengthtrace,minlengthmother,maxthreshold,minthreshold,maxpos,minpos)
%hist(max(signals,[],2),100);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
signals=tracedata(:,:,channel);
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shorttraces=tracestats(:,3)<=minlengthtrace;
shorttraces=tracestats(:,3)<minlengthtrace | motherstats(:,3)<minlengthmother;

nosignal=ones(numtraces,1)*NaN;
switch maxpos
    case 0
        nosignal=max(signals,[],2)<maxthreshold;
    case 1
        for i=1:numtraces
            nosignal(i)=signals(i,tracestats(i,1))<maxthreshold;
        end
    case 2
        for i=1:numtraces
            nosignal(i)=signals(i,tracestats(i,2))<maxthreshold;
        end
end

%strange=min(signals,[],2)>minthreshold;
strange=ones(numtraces,1)*NaN;
switch minpos
    case 0
        strange=min(signals,[],2)>minthreshold;
    case 1
        for i=1:numtraces
            strange(i)=min(signals(i,motherstats(i,1):motherstats(i,2)))>minthreshold;
        end
    case 2
        for i=1:numtraces
            strange(i)=min(signals(i,tracestats(i,1):tracestats(i,2)))>minthreshold;
        end
end

badtraces=shorttraces | nosignal | strange;
signals=tracedata(:,:,channel); %return raw signal (not smoothened)