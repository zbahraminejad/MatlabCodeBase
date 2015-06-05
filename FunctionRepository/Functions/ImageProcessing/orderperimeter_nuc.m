function order=orderperimeter_nuc(coorset)
numpoints=size(coorset,1);
order=zeros(numpoints,1);
order(1)=1;

%%% determine clockwise direction from r(1),c(1) %%%%%%%%%%%%%%%%%%
compass=[1,0;1,-1;0,-1;-1,-1;-1,0;-1,1;0,1;1,1];
buffer=zeros(size(coorset));
x=coorset(1,1); y=coorset(1,2);
buffer(:,1)=coorset(:,1)-x;
buffer(:,2)=coorset(:,2)-y;
buffer(1,:)=[10000,10000];
adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
%if numel(adjidx)~=2
%    fprintf('Other than 2 neighbors detected!\n');
%end
[~,idx]=sortrows(buffer(adjidx,:),[-2 1]);   %sort first by descending y then ascending x
adjidx=adjidx(idx(1),:);
%%% order ring coordinates contiguously %%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=2:numpoints-1
    %fprintf('pp = %0.0f\n',pp);
    order(pp)=adjidx;
    x=buffer(adjidx,1); y=buffer(adjidx,2);
    buffer(:,1)=buffer(:,1)-x;
    buffer(:,2)=buffer(:,2)-y;
    buffer(adjidx,:)= [10000,10000];
    adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
    if numel(adjidx)~=1
        %fprintf('Other than 1 neighbor detected!\n');
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
order(numpoints)=adjidx;

end