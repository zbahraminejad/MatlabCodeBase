function traces=normalizetraces(traces)
for i=1:size(traces,1)
    maxval=max(traces(i,:));
    traces(i,:)=traces(i,:)/maxval;
end
end