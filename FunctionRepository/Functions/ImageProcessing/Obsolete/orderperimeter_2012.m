function orderedset=orderperimeter(set)
orderedset=zeros(size(set));
orderedset(1,:) = set(1,:);
numpoints=size(set,1);
%%% determine clockwise direction from r(1),c(1) %%%%%%%%%%%%%%%%%%
compass=[1,0;1,-1;0,-1;-1,-1;-1,0;-1,1;0,1;1,1];
buffer=zeros(size(set));
x=set(1,1); y=set(1,2);
buffer(:,1)=set(:,1)-x;
buffer(:,2)=set(:,2)-y;
buffer(1,:)=[10000,10000];
%{
adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
if numel(adjidx)~=2
    fprintf('Other than 2 neighbors detected!\n');
end
poles=find(ismember(compass,buffer(adjidx,:),'rows'));
side=compass(round(mean(poles)),:);
nucpres = nucmask(y+side(2),x+side(1));
if nucpres
    pole = poles(1);
else
    pole = poles(2);
end
adjidx=find(buffer(:,1)==compass(pole,1) & buffer(:,2)==compass(pole,2));
%}
adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
if numel(adjidx)~=2
    fprintf('Other than 2 neighbors detected!\n');
end
[~,idx]=sortrows(buffer(adjidx,:),-2);
adjidx=adjidx(idx(1),:);

%%% order ring coordinates contiguously %%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=2:numpoints-1
    orderedset(pp,:)=set(adjidx,:);
    x=buffer(adjidx,1); y=buffer(adjidx,2);
    buffer(:,1)=buffer(:,1)-x;
    buffer(:,2)=buffer(:,2)-y;
    buffer(adjidx,:)= [10000,10000];
    adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
    if numel(adjidx)~=1
        pole=find(compass(:,1)==x & compass(:,2)==y);
        antipole=pole+4;
        priority=[antipole+1:antipole+7]';
        priority=mod(priority,8); priority(priority==0)=8;
        priority=compass(priority,:);
        poles=find(ismember(priority,buffer(adjidx,:),'rows'));
        x=priority(poles(1),1); y=priority(poles(1),2);
        adjidx=find(buffer(:,1)==x & buffer(:,2)==y);
    end
end
orderedset(numpoints,:)=set(adjidx,:);

end