function [x,y] = calculatejitters_cellmatch(wellsss,badframes,height,width)
totalframes = size(wellsss,3);
frameidx = find(badframes==1);
EF = totalframes-length(frameidx);
midframe=round(totalframes/2); midframe=midframe+ismember(midframe,frameidx);
sizediffiqr=iqr(wellsss{midframe}(:,4));
x=zeros(1,EF); y=zeros(1,EF);
winrad=100;
g = 1;
for f=1:EF
    g=g+ismember(f,frameidx);
    xx=wellsss{g}(:,1);
    yy=wellsss{g}(:,2);
    leavein=find(xx>winrad & xx<width-winrad & yy>winrad & yy<height-winrad);
    xx=xx(leavein);
    yy=yy(leavein);
    cellsize=wellsss{g}(leavein,4);
    cellnum=numel(xx);
    sizediff1=zeros(cellnum,1);
    sizediff2=zeros(cellnum,1);
    match=zeros(cellnum,1);
    xdiff=zeros(cellnum,1);
    ydiff=zeros(cellnum,1);
    if f>1
        for i=1:cellnum
            neighbors=find(abs(wellsss{gprev}(:,1)-xx(i))<winrad & abs(wellsss{gprev}(:,2)-yy(i))<winrad);
            sizediff=abs(wellsss{gprev}(neighbors,4)-cellsize(i))/sizediffiqr;
            if numel(neighbors)==1
                sizediff1(i)=sizediff;
                match(i)=sizediff<0.2;
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
            lh2=sizediff(nextbest)-sizediff(samecell)>0.2;
            match(i)=lh1 && lh2;
            xdiff(i)=wellsss{gprev}(neighbors(samecell),1)-xx(i);
            ydiff(i)=wellsss{gprev}(neighbors(samecell),2)-yy(i);
        end
        xdiff=xdiff(match==1);
        ydiff=ydiff(match==1);
        [xmu,xs]=normfit(xdiff);
        [ymu,ys]=normfit(ydiff);
        keep=find(xdiff>xmu-xs & xdiff<xmu+xs & ydiff>ymu-ys & ydiff<ymu+ys);
        xdiff=xdiff(keep);
        ydiff=ydiff(keep);
        x(f)=x(f-1)+mean(xdiff);
        y(f)=y(f-1)+mean(ydiff);
        %fprintf('totalcells = %0.0f, goodcells = %0.0f, x = %0.2f, y = %0.2f\n',cellnum,numel(xdiff),x(f),y(f));
    end
    gprev=g;
    g = g+1;
end