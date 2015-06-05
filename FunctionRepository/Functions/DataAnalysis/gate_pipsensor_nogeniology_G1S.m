function badtraces=gate_pipsensor(data,tracestats,negthreshold,maxthreshold,minthreshold,initialwindow)
%hist(max(signals,[],2),100);
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(data,1);
signals=data;
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
negmin=min(signals,[],2)<negthreshold;

toolow=zeros(numtraces,1);

for i=1:numtraces
    if tracestats(i,3)>initialwindow
        toolow(i)=max(signals(i,tracestats(i,1):tracestats(i,1)+initialwindow-1))<maxthreshold;
    end
end

toohigh=ones(numtraces,1)*NaN;
for i=1:numtraces
	toohigh(i)=min(signals(i,tracestats(i,1):tracestats(i,2)))>minthreshold;
end

badtraces=negmin | toolow | toohigh;