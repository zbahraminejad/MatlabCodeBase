function [orderedidx,sortedPOI]=GetPointsOfInterest(signal,POItraces,stats,minlength,alignoption)
numtraces=size(POItraces,1);
sampleid=(1:numtraces)';
%%% detect onset of CDK2 activity %%%%%%%%%%%
relative=1; %returns Cdk2start time relative to mitosis
switch signal
    case 'Cdk2'
        [~,POI,~,badtraces]=getCdk2features(POItraces,stats,minlength,relative);
    case 'Geminin'
        [POI,badtraces]=getGemininFeatures_test2(POItraces,stats,minlength,relative);
end
sampleid(badtraces>0)=[];
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if alignoption==1
    [~,idx]=sort(POI(sampleid));
elseif alignoption==2
    [~,idx]=sort(stats(sampleid,1));
end
orderedidx=sampleid(idx);
sortedPOI=POI(orderedidx);
% realidx=sampleid(~isnan(POI(sampleid)));
% realPOI=POI(realidx);