function [xdiff,ydiff] = getjitter_areamass(centers,areas,masses,prevdata,dims)
height=dims(1); width=dims(2);
winrad=100;

xcur=centers(:,1); ycur=centers(:,2);
xprev=prevdata(:,1); yprev=prevdata(:,2); areasprev=prevdata(:,3); massesprev=prevdata(:,4);
leavein=find(xcur>winrad & xcur<width-winrad & ycur>winrad & ycur<height-winrad);
xcur=xcur(leavein);
ycur=ycur(leavein);
numcells=numel(xcur);
areadiff1=zeros(numcells,1);
areadiff2=zeros(numcells,1);
massdiff1=zeros(numcells,1);
massdiff2=zeros(numcells,1);
match=zeros(numcells,1);
xdiff=zeros(numcells,1);
ydiff=zeros(numcells,1);
for i=1:numcells
    neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    areadiff=abs(areasprev(neighbors)-areas(i))/areas(i);
    massdiff=abs(massesprev(neighbors)-masses(i))/masses(i);
    if isempty(neighbors)
    elseif numel(neighbors)==1
        areadiff1(i)=areadiff;
        massdiff1(i)=massdiff;
        match(i)=areadiff<0.05 & massdiff<0.05;
        xdiff(i)=xprev(neighbors)-xcur(i);
        ydiff(i)=yprev(neighbors)-yprev(i);
    else
        [~,closestidx]=sort(2*massdiff+areadiff);
        samecell=closestidx(1);
        nextbest=closestidx(2);
        areadiff1(i)=areadiff(samecell);
        areadiff2(i)=areadiff(nextbest);
        massdiff1(i)=massdiff(samecell);
        massdiff2(i)=massdiff(nextbest);
        lh1=areadiff1(i)<0.05 && massdiff1(i)<0.05;
        %lh2=(areadiff2(i)-areadiff1(i))>0.15 || (massdiff2(i)-massdiff1(i))>0.15;
        match(i)=lh1;
        xdiff(i)=xprev(neighbors(samecell))-xcur(i);
        ydiff(i)=yprev(neighbors(samecell))-ycur(i);
    end
end
xdiff=xdiff(match==1);
ydiff=ydiff(match==1);
[xmu,xs]=normfit(xdiff);
[ymu,ys]=normfit(ydiff);
keep=find(xdiff>xmu-xs & xdiff<xmu+xs & ydiff>ymu-ys & ydiff<ymu+ys);
xdiff=mean(xdiff(keep));
ydiff=mean(ydiff(keep));
%fprintf('totalcells = %0.0f, goodcells = %0.0f, x = %0.2f, y = %0.2f\n',cellnum,numel(xdiff),x(f),y(f));