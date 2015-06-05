function [Cdk2inc,Cdk2low]=categorizeCdk2fate(tracesCdk2,daughterstats,minlengthdaughter)
earlytime=15; %default 15
numtraces=size(tracesCdk2,1);
earlyval=ones(numtraces,1)*NaN;
lateval=ones(numtraces,1)*NaN;
maxval=ones(numtraces,1)*NaN;
minval=ones(numtraces,1)*NaN;
for i=1:numtraces
    earlyval(i)=tracesCdk2(i,daughterstats(i,1)+earlytime-2:daughterstats(i,1)+earlytime-1); 
    lateval(i)=tracesCdk2(i,daughterstats(i,1)+minlengthdaughter-1);
    maxval(i)=max(tracesCdk2(i,daughterstats(i,1)+earlytime:daughterstats(i,2))); %max from minlength-end
    minval(i)=min(tracesCdk2(i,daughterstats(i,1)+earlytime:daughterstats(i,2)));
end
%figure,hist(earlyval,0.3:0.02:1.5);
earlycutoff=0.6; %default 0.55
latecutoff=earlycutoff+0.01*(minlengthdaughter-earlytime);
Cdk2inc=earlyval>earlycutoff & lateval>latecutoff & minval>earlycutoff;
Cdk2low=earlyval<earlycutoff & maxval<earlycutoff;
end