function order=orderperimeter_fast(coorset)
numpoints=size(coorset,1);
order=zeros(numpoints,1);
order(1)=1;
xorg=coorset(1,1); yorg=coorset(1,2);
coorset(1,:)=[0 0];
idx=find(coorset(:,1)==xorg+1 & coorset(:,2)==yorg-1);
curpoint=coorset(idx,:);
coorset(idx,:)=[0 0];
order(2)=idx;
lastidx=find(ismember(coorset(:,1),[xorg xorg+1]) & coorset(:,2)==yorg+1);
c=3;
while idx~=lastidx
    idx=find(ismember(coorset(:,1),curpoint(1)-1:curpoint(1)+1) & ismember(coorset(:,2),curpoint(2)-1:curpoint(2)+1),1,'first');
    curpoint=coorset(idx,:);
    coorset(idx,:)=[0 0];
    order(c)=idx;
    c=c+1;
end
end