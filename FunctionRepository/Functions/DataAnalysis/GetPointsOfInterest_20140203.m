function [datasorted,POI]=GetPointsOfInterest(signal,POItraces,datatotal,stats,minlength,alignoption)
numtraces=size(POItraces,1);
sampleid=(1:numtraces)';
%%% detect onset of CDK2 activity %%%%%%%%%%%
relative=1; %returns Cdk2start time relative to mitosis
switch signal
    case 'Cdk2'
        [~,POI,~,badtraces]=getCdk2features(POItraces,stats,minlength,relative);
    case 'Geminin'
        [POI,badtraces]=getGemininFeatures(POItraces,stats,minlength,relative);
end
POI(isnan(POI))=0;
sampleid(badtraces>0)=[];
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alignoption==1
    [~,idx]=sort(POI(sampleid));
elseif alignoption==2
    [~,idx]=sort(stats(sampleid,1));
end
ordered=sampleid(idx);
POI=POI(ordered);
datasorted=datatotal(ordered,:);