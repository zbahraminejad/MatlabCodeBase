function nans=gate_nan(tracedata,tracestats)
%hist(min(signals,[],2),0:10:1001); xlim([0 1000]);
%hist(max(signals,[],2),0:10:2010); xlim([0 2000]);
%hist(log2(min(signals,[],2)),0:0.1:15.1); xlim([0 15]);

%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
nans=ones(numtraces,1)*NaN;
for i=1:numtraces
    nans(i)=sum(isnan(tracedata(i,tracestats(i,1):tracestats(i,2),1)));
end
nans=nans>0;