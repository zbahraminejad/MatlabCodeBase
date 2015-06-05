function [time,datatotal]=alignmotherdaughter2(alltraces,allstats,allmotherstats,alignoption,historyoption,motherperc,daughterperc,framesperhr)
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
elseif alignoption==3
    numtraces=size(alltraces1,1);
    IFalignedtraces1=ones(numtraces,numframes)*NaN;
    IFalignedtraces2=ones(numtraces,numframes)*NaN;
    for i=1:numtraces
        IFalignedtraces1(i,numframes-allstats(i,3)+1:numframes)=alltraces1(i,allstats(i,1):allstats(i,2));
        IFalignedtraces2(i,numframes-allstats(i,3)+1:numframes)=alltraces2(i,allstats(i,1):allstats(i,2));
    end
    maxlength=prctile(allstats(:,3),90);
    IFalignedtraces1=IFalignedtraces1(:,numframes-maxlength+1:numframes);
    IFalignedtraces2=IFalignedtraces2(:,numframes-maxlength+1:numframes);
end