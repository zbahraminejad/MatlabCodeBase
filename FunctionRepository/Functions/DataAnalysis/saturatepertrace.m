function saturatetraces=saturatepertrace(traces)
[numtraces,numframes]=size(traces);
saturatetraces=ones(numtraces,numframes);
for i=1:numtraces
    tracenan=isnan(traces(i,:));
    tracesat=imadjust(mat2gray(traces(i,:)));
    tracesat(tracenan)=NaN;
    saturatetraces(i,:)=tracesat;
end