function traces=normalizetraces_3(traces,tracestats,normoption)
smoothtraces=traces;
minval=ones(size(traces,1),1);
maxval=ones(size(traces,1),1);
for i=1:size(traces,1)
    realframes=find(~isnan(traces(i,:)));
    smoothtraces(i,realframes)=smooth(traces(i,realframes));
    minval(i)=min(smoothtraces(i,:));
    maxval(i)=max(smoothtraces(i,:));
end
medianmax=nanmedian(maxval);
max90=prctile(maxval,90);
for i=1:size(traces,1)
    if normoption==0
        traces(i,:)=traces(i,:)/maxval(i);
    elseif normoption==1 && minval(i)>0.5*maxval(i)
        traces(i,:)=traces(i,:)/medianmax;
    elseif normoption==1 && minval(i)<=0.5*maxval(i)
        traces(i,:)=traces(i,:)/maxval(i);
    elseif normoption==2
        traces(i,:)=traces(i,:)/smoothtraces(i,tracestats(i,1));
    elseif normoption==3
        traces(i,:)=traces(i,:)/max90;
    end
end
end