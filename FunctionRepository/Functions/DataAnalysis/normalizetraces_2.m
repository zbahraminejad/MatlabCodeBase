function traces=normalizetraces_1(traces,rangeoption)
minval=ones(size(traces,1),1);
maxval=ones(size(traces,1),1);
for i=1:size(traces,1)
    realframes=find(~isnan(traces(i,:)));
    traces(i,realframes)=smooth(traces(i,realframes));
    minval(i)=min(traces(i,:));
    maxval(i)=max(traces(i,:));
end
medianmax=nanmedian(maxval);
for i=1:size(traces,1)
    if rangeoption && minval(i)>0.9*maxval(i)
        traces(i,:)=traces(i,:)/medianmax;
    else
        traces(i,:)=traces(i,:)/maxval(i);
    end
end
end