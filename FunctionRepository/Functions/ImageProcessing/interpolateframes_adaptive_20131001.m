function [tracedata,jitters]=interpolateframes_adaptive(tracedata,jitters,badframes)
if badframes>0
    totalframes=size(tracedata,1);
    for i=badframes'
        if i==1
            tracedata{1}=tracedata{2};
            jitters(1,:)=jitters(2,:);
        elseif i==totalframes
            tracedata{i}=tracedata{i-1};
            jitters(i,:)=jitters(i-1,:);
        else
            prevnum=find(~isnan(tracedata{i-1}(:,1)),1,'last');
            tracedata{i}=tracedata{i-1}(1:prevnum,:)+tracedata{i+1}(1:prevnum,:)/2;
            jitters(i,:)=mean([jitters(i-1,:);jitters(i+1,:)]);
        end
    end
end

end