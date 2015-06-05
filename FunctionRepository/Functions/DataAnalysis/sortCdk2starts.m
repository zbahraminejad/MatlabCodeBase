function [datasorted,Cdk2start]=sortCdk2starts(alltracesCdk2,datatotal,alldaughterstats,minlengthdaughter,alignoption)
numtraces=size(alltracesCdk2,1);
sampleid=(1:numtraces)';
%%% detect onset of CDK2 activity %%%%%%%%%%%
relative=1; %returns Cdk2start time relative to mitosis
[~,Cdk2start,~,badtraces]=getCdk2features(alltracesCdk2,alldaughterstats,minlengthdaughter,relative);
Cdk2start(isnan(Cdk2start))=0;
sampleid(badtraces>0)=[];
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alignoption==1
    [~,idx]=sort(Cdk2start(sampleid));
elseif alignoption==2
    [~,idx]=sort(alldaughterstats(sampleid,1));
end
ordered=sampleid(idx);
Cdk2start=Cdk2start(ordered);
datasorted=datatotal(ordered,:);