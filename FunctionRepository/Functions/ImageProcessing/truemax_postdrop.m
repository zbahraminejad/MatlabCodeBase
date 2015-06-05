function boundarycheck = truemax_postdrop(actualtime,prevpoint,nextpoint,smoothcdt1)
lowercheck = sphasetime<=curpoint+2;
uppercheck = sphasetime>=nextpoint-20;
if uppercheck==0
    nextmin = min(smoothcdt1(sphasetime:nextpoint));
    if nextmin>smoothcdt1(sphasetime)*0.5 && smoothangie(sphasetime)<0.75
        uppercheck=1;
    end
end
boundarycheck = lowercheck || uppercheck;
end
