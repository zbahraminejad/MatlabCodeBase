function badtraces=gate_Geminin_9_mother(data,tracestats,maxthreshold,minthreshold,initialwindow)
%hist(max(signals,[],2),100);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(data,1);
signals=data;
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toolow=ones(numtraces,1)*NaN;
for i=1:numtraces
    toolow(i)=max(signals(i,tracestats(i,1):tracestats(i,2)))<maxthreshold;
end

toohigh=zeros(numtraces,1);
for i=1:numtraces
    if tracestats(i,3)>initialwindow
        toohigh(i)=min(signals(i,tracestats(i,1):tracestats(i,1)+initialwindow-1))>minthreshold;
    end
end

badtraces=toolow | toohigh;