function [tracedata,jitters]=interpolateframes_MZ(tracedata,jitters,badframe)

prevnum=find(~isnan(tracedata(:,badframe-1,1)),1,'last');
tracedata(1:prevnum,badframe,:)=(tracedata(1:prevnum,badframe-1,:)+tracedata(1:prevnum,badframe+1,:))/2;
jitters(badframe,:)=mean([jitters(badframe-1,:);jitters(badframe+1,:)]);



end