function boundarycheck = truemax(actualtime,prevpoint,nextpoint)
lowercheck = actualtime<=prevpoint+2;   %check left side
uppercheck = actualtime>=nextpoint-20;   %check right side
boundarycheck = lowercheck || uppercheck;
end
