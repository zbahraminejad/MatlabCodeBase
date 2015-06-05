function [Cdk2inc,Cdk2low]=categorizeCdk2fate_3(tracesCdk2,daughterstats,minlengthdaughter)
earlytime=14; %default 14
numtraces=size(tracesCdk2,1);
earlyval=ones(numtraces,1)*NaN;
lateval=ones(numtraces,1)*NaN;
maxval=ones(numtraces,1)*NaN;
meanval=ones(numtraces,1)*NaN;
for i=1:numtraces
    earlyval(i)=tracesCdk2(i,daughterstats(i,1)+earlytime); 
    lateval(i)=tracesCdk2(i,daughterstats(i,1)+minlengthdaughter-1);
    maxval(i)=max(tracesCdk2(i,daughterstats(i,1)+earlytime:daughterstats(i,2))); %max from minlength-end
    %minval(i)=min(tracesCdk2(i,daughterstats(i,1)+earlytime:daughterstats(i,2)));
    meanval(i)=mean(tracesCdk2(i,daughterstats(i,1)+earlytime:daughterstats(i,2)));
end
%figure,hist(earlyval,0:0.05:1.6); xlim([0.1 1.5]);
earlycutoff=0.55; %default 0.55
latecutoff=earlycutoff+0.01*(minlengthdaughter-earlytime);
Cdk2inc=earlyval>earlycutoff & lateval>latecutoff & meanval>earlycutoff;
Cdk2low=earlyval<earlycutoff & maxval<earlycutoff;
end