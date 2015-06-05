function traces=normalizetraces(traces,stats)
for i=1:size(traces,1)
    maxval=max(traces(i,stats(i,1):stats(i,2)));
    traces(i,:)=traces(i,:)/maxval;
end
end