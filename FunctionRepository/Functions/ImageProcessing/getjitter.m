function [xdiff,ydiff] = getjitter(centers,centersprev,areas,areasprev,dims)
height=dims(1); width=dims(2);
sizediffiqr=iqr(areas);
winrad=100;

xcur=centers(:,1); ycur=centers(:,2);
xprev=centersprev(:,1); yprev=centersprev(:,2);
leavein=find(xcur>winrad & xcur<width-winrad & ycur>winrad & ycur<height-winrad);
xcur=xcur(leavein);
ycur=ycur(leavein);
numcells=numel(xcur);
sizediff1=zeros(numcells,1);
sizediff2=zeros(numcells,1);
match=zeros(numcells,1);
xdiff=zeros(numcells,1);
ydiff=zeros(numcells,1);
for i=1:numcells
    neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    sizediff=abs(areasprev(neighbors)-areas(i))/sizediffiqr;
    if numel(neighbors)==1
        sizediff1(i)=sizediff;
        match(i)=sizediff<0.2;
        xdiff(i)=xprev(neighbors)-xcur(i);
        ydiff(i)=yprev(neighbors)-ycur(i);
        continue;
    end
    if isempty(neighbors)
        continue;
    end
    [~,closestidx]=sort(sizediff);
    samecell=closestidx(1);
    nextbest=closestidx(2);
    sizediff1(i)=sizediff(samecell);
    sizediff2(i)=sizediff(nextbest);
    lh1=sizediff(samecell)<0.2;
    lh2=(sizediff(nextbest)-sizediff(samecell))>0.2;
    match(i)=lh1 && lh2;
    xdiff(i)=xprev(neighbors(samecell))-xcur(i);
    ydiff(i)=yprev(neighbors(samecell))-ycur(i);
end
xdiff=xdiff(match==1);
ydiff=ydiff(match==1);
[xmu,xs]=normfit(xdiff);
[ymu,ys]=normfit(ydiff);
keep=find(xdiff>xmu-xs & xdiff<xmu+xs & ydiff>ymu-ys & ydiff<ymu+ys);
xdiff=mean(xdiff(keep));
ydiff=mean(ydiff(keep));
%fprintf('totalcells = %0.0f, goodcells = %0.0f, x = %0.2f, y = %0.2f\n',cellnum,numel(xdiff),x(f),y(f));