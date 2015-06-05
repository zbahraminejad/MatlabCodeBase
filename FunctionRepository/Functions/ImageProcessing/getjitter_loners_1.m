function [jitx,jity,singlecount] = getjitter_loners(centers,prevdata,winrad,dims)
height=dims(1); width=dims(2);

xcur=centers(:,1); ycur=centers(:,2);
xprev=prevdata(:,1); yprev=prevdata(:,2);
leavein=find(xcur>winrad & xcur<width-winrad & ycur>winrad & ycur<height-winrad);
xcur=xcur(leavein);
ycur=ycur(leavein);
numcells=numel(xcur);
xdiff=ones(numcells,1)*NaN;
ydiff=ones(numcells,1)*NaN;
singlecount=0;
for i=1:numcells
    neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    if numel(neighbors)==1
        singlecount=singlecount+1;
        xdiff(i)=xprev(neighbors)-xcur(i);
        ydiff(i)=yprev(neighbors)-ycur(i);
    end
end
jitx=mean(xdiff(~isnan(xdiff)));
jity=mean(ydiff(~isnan(ydiff)));
%[xmu,xs]=normfit(xdiff);
%[ymu,ys]=normfit(ydiff);
% keep=find(xdiff>xmu-xs & xdiff<xmu+xs & ydiff>ymu-ys & ydiff<ymu+ys);
% xdiff=xdiff(keep);
% ydiff=ydiff(keep);

%fprintf('totalcells = %0.0f, goodcells = %0.0f, x = %0.2f, y = %0.2f\n',cellnum,numel(xdiff),x(f),y(f));