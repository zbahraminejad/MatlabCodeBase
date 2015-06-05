function [POI,badtraces]=GetPointsOfInterest_1(signal,POItraces,stats,minlength)
switch signal
    case 'Cdk2'
        [POI,badtraces]=getCdk2features_1(POItraces,stats,minlength);
    case 'Geminin'
        [POI,badtraces]=getGemininfeatures_1(POItraces,stats);
end