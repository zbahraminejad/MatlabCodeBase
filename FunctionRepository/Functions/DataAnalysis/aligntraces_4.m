function [time,datatotal]=aligntraces_4(alltraces,alignPOI,allstats,allmotherstats,daughteroption)
[numtraces,numframes]=size(alltraces);
historyoption=1; %1:motherdata 2:all ancestry
if daughteroption==2
    historyoption=2; %must do this if no mother data exists
end
tracedataleft=ones(numtraces,numframes)*NaN;
tracedataright=ones(numtraces,numframes)*NaN;
lengthleft=ones(numtraces,1)*NaN; lengthright=lengthleft;
for i=1:numtraces
    if historyoption==1 %only show mother data
        firstframe=allmotherstats(i,1);
    elseif historyoption==2 %show entire ancestry
        firstframe=find(~isnan(alltraces(i,:)),1,'first');
    end
    if isnan(alignPOI(i))
        alignPOI(i)=allstats(i,2); %NaN values set to end of trace only for aligning
    end
    lengthleft(i)=alignPOI(i)-firstframe+1;
    tracedataleft(i,numframes-lengthleft(i)+1:numframes)=alltraces(i,firstframe:alignPOI(i));
    lengthright(i)=allstats(i,2)-alignPOI(i)+1;
    tracedataright(i,1:lengthright(i))=alltraces(i,alignPOI(i):allstats(i,2));
end
%lengthleftmax=floor(prctile(allmotherstats(:,3),90));
lengthleftmax=floor(prctile(lengthleft,90));
tracedataleft=tracedataleft(:,numframes-lengthleftmax+1:numframes);
%rightlengthmax=floor(prctile(allstats(:,3),90));
rightlengthmax=floor(prctile(lengthright,90));
tracedataright=tracedataright(:,1:rightlengthmax);
datatotal=[tracedataleft tracedataright];
time=(-lengthleftmax+1:rightlengthmax);