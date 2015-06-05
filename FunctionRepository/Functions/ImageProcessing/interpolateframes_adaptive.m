function [tracedata,jitters]=interpolateframes_adaptive(tracedata,jitters,badframes)

totalframes=size(tracedata,1);
badframeid=find(badframes);
for i=badframeid'
    if i==1
        firstgood=find(badframes==0,1,'first');
        tracedata{i}=tracedata{firstgood};
        jitters(i,:)=jitters(firstgood,:);
    elseif i==totalframes || badframes(i+1)==1
        tracedata{i}=tracedata{i-1};
        jitters(i,:)=jitters(i-1,:);
    else
        prevnum=find(~isnan(tracedata{i-1}(:,1)),1,'last');
        tracedata{i}=tracedata{i-1}(1:prevnum,:)+tracedata{i+1}(1:prevnum,:)/2;
        jitters(i,:)=mean([jitters(i-1,:);jitters(i+1,:)]);
    end
    badframes(i)=0;
end

end