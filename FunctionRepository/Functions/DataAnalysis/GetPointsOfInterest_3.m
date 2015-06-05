function [POI,badtraces]=GetPointsOfInterest_2(signal,traces,stats,minlength)
switch signal
    case 'Cdk2'
        %[POI,badtraces]=getCdk2features_1(traces,stats,minlength);
        %[POI,badtraces]=getCdk2features_steve(traces,stats,minlength);
        [POI,badtraces]=getCdk2features_2(traces,stats,minlength);
    case 'Geminin'
        [POI,badtraces]=getGemininfeatures_1(traces,stats);
    case 'Geminin_SR'
        [POI,badtraces]=getGemininfeatures_serumrelease(traces,stats);
    case 'PipdegFall'
        %[POI,badtraces]=getPipdegFall_NormByFirst(traces,stats,minlength);
        %[POI,badtraces]=getPipdegFall_NormByMother(traces,stats,minlength);
        %[POI,badtraces]=getPipdegFall_Dropstart(traces,stats,minlength);
        %[POI,badtraces]=getPipdegFall_Dropstart_1(traces,stats);
        %[POI,badtraces]=getPipdegFall_Dropstart_2(traces,stats,minlength);
        [POI,badtraces]=getPipdegFall_Dropstart_3(traces,stats,minlength);
    case 'PipdegRise'
        %[POI,badtraces]=getPipdegRise(traces,stats);
        [POI,badtraces]=getPipdegRise_2(traces,stats,minlength);
end