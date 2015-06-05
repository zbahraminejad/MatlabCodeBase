function normtraces=normalizepertrace(traces)
[numtraces,numframes]=size(traces);
normtraces=ones(numtraces,numframes);
for i=1:numtraces
    normtraces(i,:)=traces(i,:)/max(traces(i,:));
end