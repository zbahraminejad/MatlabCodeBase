function traces=normalizetraces_1(traces,rangeoption)
minval=min(traces,[],2);
maxval=max(traces,[],2);
medianmax=median(maxval);
for i=1:size(traces,1)
    realframes=find(~isnan(traces(i,:)));
    traces(i,realframes)=smooth(traces(i,realframes));
    if rangeoption && minval(i)>0.9*maxval(i)
        traces(i,:)=traces(i,:)/medianmax;
    else
        traces(i,:)=traces(i,:)/maxval;
    end
end
end