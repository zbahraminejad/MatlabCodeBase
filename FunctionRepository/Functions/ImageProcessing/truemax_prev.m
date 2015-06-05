function boundarycheck = truemax(cdt1,actualtime,numframes,edge,side)
lowertimebound = actualtime-20;
uppertimebound = actualtime+20;
%%% check left side %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataedge = lowertimebound<=0;
if dataedge
    lowercheck = 1;
else
    if actualtime<=3
        stillrising = 1;
    elseif strcmp(side,'left')
        stillrising = actualtime==edge && cdt1(actualtime-3)>cdt1(actualtime);
    else
        stillrising = 0;
    end
    slowfall = cdt1(actualtime)-cdt1(lowertimebound)<0.3;
    traceedge = cdt1(actualtime-5)==-10000;
    lowercheck = dataedge || stillrising || slowfall || traceedge;
end
%%% check right side %%%%%%%%%%%%%%%%%%%%%%%%%%%
dataedge = uppertimebound>numframes;
if dataedge
    uppercheck = 1;
else
    if actualtime>=numframes-2
        stillrising = 1;
    elseif strcmp(side,'right')
        stillrising = actualtime==edge && cdt1(actualtime+3)>cdt1(actualtime);
    else
        stillrising = 0;
    end
    slowfall = cdt1(actualtime)-cdt1(uppertimebound)<0.3;
    traceedge = cdt1(actualtime+5)==-10000;
    uppercheck = dataedge || stillrising || slowfall || traceedge;
end
boundarycheck = lowercheck || uppercheck;
end
