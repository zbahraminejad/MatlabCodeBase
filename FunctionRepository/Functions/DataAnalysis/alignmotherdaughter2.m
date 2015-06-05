function [time,datatotal]=alignmotherdaughter2(alltraces,allstats,allmotherstats,alignoption,drugspike,historyoption,motherperc,daughterperc,framesperhr)
[numtraces,numframes]=size(alltraces);
if alignoption==1
    datadaughter=ones(numtraces,numframes)*NaN;
    datamother=ones(numtraces,numframes)*NaN;
    for i=1:numtraces
        datadaughter(i,1:allstats(i,3))=alltraces(i,allstats(i,1):allstats(i,2));
        lastframe=allmotherstats(i,2);
        if historyoption==1 %only show mother data
            firstframe=allmotherstats(i,1);
        elseif historyoption==2 %show entire ancestry
            firstframe=find(~isnan(alltraces(i,:)),1,'first');
        end
        duration=lastframe-firstframe+1;
        datamother(i,numframes-duration+1:numframes)=alltraces(i,firstframe:lastframe);
    end
    daughtermax=floor(prctile(allstats(:,3),daughterperc));
    datadaughter=datadaughter(:,1:daughtermax);
    mothermax=floor(prctile(allmotherstats(:,3),motherperc));
    datamother=datamother(:,numframes-mothermax+1:numframes);
    datatotal=[datamother datadaughter];
    time=(-mothermax+1:daughtermax)/framesperhr;
elseif alignoption==2
    daughtermax=prctile(allstats(:,2),daughterperc);
    datatotal=ones(numtraces,daughtermax)*NaN;
    for i=1:numtraces
        if historyoption==1
            datatotal(i,allmotherstats(i,1):daughtermax)=alltraces(i,allmotherstats(i,1):daughtermax);
        elseif historyoption==2
            datatotal(i,1:daughtermax)=alltraces(i,daughtermax);
        end
    end
    time=(1-drugspike:daughtermax-drugspike)/framesperhr;
end