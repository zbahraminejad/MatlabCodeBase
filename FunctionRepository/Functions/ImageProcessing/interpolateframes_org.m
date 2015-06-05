function [bestsp,best_rc,x,y]=interpolateframes(bestsp,best_rc,x,y,badframes,fmax)

if badframes>0
    tempsp=cell(1,1,fmax); tempx=zeros(1,fmax); tempy=zeros(1,fmax);
    j=0; prev=1;
    for i=badframes
        if i==1
            tempsp{1}=bestsp{1};
            tempx(1)=x(1); tempy(1)=y(1);
        elseif i==fmax
            tempsp(prev:i-1)=bestsp(prev-j:i-j-1);
            tempsp{fmax}=bestsp{end};
            tempx(prev:i-1)=x(prev-j:i-j-1);
            tempx(fmax)=x(end);
            tempy(prev:i-1)=y(prev-j:i-j-1);
            tempy(fmax)=y(end);
        else
            tempsp(prev:i-1)=bestsp(prev-j:i-j-1);
            tempsp{i}=(bestsp{i-j-1}+bestsp{i-j}(1:size(bestsp{i-j-1},1),:))/2;    %average of each value for preceding and following frames
            tempx(prev:i-1)=x(prev-j:i-j-1);
            tempx(i)=(x(i-j-1)+x(i-j))/2;
            tempy(prev:i-1)=y(prev-j:i-j-1);
            tempy(i)=(y(i-j-1)+y(i-j))/2;
        end
        startidx=find(best_rc(:,1)>=i);
        endidx=find(best_rc(:,3)>=i);
        best_rc(startidx,1)=best_rc(startidx,1)+1;
        best_rc(endidx,3)=best_rc(endidx,3)+1;
        
        j=j+1;
        prev=i+1;
    end
    if prev<=fmax
        tempsp(prev:fmax)=bestsp(prev-j:end);
        tempx(prev:fmax)=x(prev-j:end);
        tempy(prev:fmax)=y(prev-j:end);
    end
    bestsp=tempsp;
    x=tempx;
    y=tempy;
end

end