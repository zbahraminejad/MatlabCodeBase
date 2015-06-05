function [time,datatotal]=alignmotherdaughter(alltraces,alldaughterstats,allmotherstats,daughtermax,mothermax,framesperhr)
[numtraces,numframes]=size(alltraces);
datadaughter=ones(numtraces,numframes)*NaN;
datamother=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    datadaughter(i,1:alldaughterstats(i,3))=alltraces(i,alldaughterstats(i,1):alldaughterstats(i,2));
    datamother(i,numframes-allmotherstats(i,3)+1:numframes)=alltraces(i,allmotherstats(i,1):allmotherstats(i,2));
end
datadaughter=datadaughter(:,1:daughtermax);
datamother=datamother(:,numframes-mothermax+1:numframes);
datatotal=[datamother datadaughter];
time=(-mothermax+1:daughtermax)/framesperhr;