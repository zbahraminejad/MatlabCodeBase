function [xdiff,ydiff] = getjitter_mass(centers,masses,prevdata,dims)
height=dims(1); width=dims(2);
winrad=100;

xcur=centers(:,1); ycur=centers(:,2);
xprev=prevdata(:,1); yprev=prevdata(:,2); massesprev=prevdata(:,3);
leavein=find(xcur>winrad & xcur<width-winrad & ycur>winrad & ycur<height-winrad);
xcur=xcur(leavein);
ycur=ycur(leavein);
numcells=numel(xcur);
massdiff1=zeros(numcells,1);
massdiff2=zeros(numcells,1);
match=zeros(numcells,1);
xdiff=zeros(numcells,1);
ydiff=zeros(numcells,1);
for i=1:numcells
    neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    massdiff=abs(massesprev(neighbors)-masses(i))/masses(i);
    if numel(neighbors)==1
        massdiff1(i)=massdiff;
        match(i)=massdiff<0.1;
        continue;
    end
    if isempty(neighbors)
        continue;
    end
    [~,closestidx]=sort(massdiff);
    samecell=closestidx(1);
    nextbest=closestidx(2);
    massdiff1(i)=massdiff(samecell);
    massdiff2(i)=massdiff(nextbest);
    lh1=massdiff(samecell)<0.1;
    %lh2=massdiff(nextbest)-massdiff(samecell)>0.2;
    lh2=massdiff(nextbest)>0.2;
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