function [time,datatotal]=aligntraces_2(alltraces,alignPOI,allstats,allmotherstats)
[numtraces,numframes]=size(alltraces);
historyoption=1; %1:motherdata 2:all ancestry
tracedataleft=ones(numtraces,numframes)*NaN;
tracedataright=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    if historyoption==1 %only show mother data
        firstframe=allmotherstats(i,1);
    elseif historyoption==2 %show entire ancestry
        firstframe=find(~isnan(alltraces(i,:)),1,'first');
    end
    lengthleft=alignPOI(i)-firstframe+1;
    tracedataleft(i,numframes-lengthleft+1:numframes)=alltraces(i,firstframe:alignPOI(i));
    lengthright=allstats(i,2)-alignPOI(i)+1;
    tracedataright(i,1:lengthright)=alltraces(i,alignPOI(i):allstats(i,2));
end
lengthleftmax=floor(prctile(allmotherstats(:,3),90));
tracedataleft=tracedataleft(:,numframes-lengthleftmax+1:numframes);
rightlengthmax=floor(prctile(allstats(:,3),90));
tracedataright=tracedataright(:,1:rightlengthmax);
datatotal=[tracedataleft tracedataright];
time=(-lengthleftmax+1:rightlengthmax);