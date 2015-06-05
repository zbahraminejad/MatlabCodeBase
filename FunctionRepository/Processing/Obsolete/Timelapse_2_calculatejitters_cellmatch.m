%function Timelapse_2_calculatejitters(row,col,site)
%row=3;col=3;site=1;
row='D';col='05';site='4';
codepath='H:\Documents\MATLAB\Imaging\';
cd([codepath,'Functions\']); %directory for all subfunctions

projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath = '2013-06-07_p21_cy2_deletions\Experiment_20130715\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
datadir = ([projectpath,experimentpath,'Data_Test\']);
rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
maskdir = [imagepath,experimentpath,'Mask_Test\',shot,'\'];

timejittercalc=tic;
load([datadir,'wellsss_',shot,'.mat'],'wellsss','badframes');
totalframes = size(wellsss,3);
frameidx = find(badframes==1);
EF = totalframes-length(frameidx);
tempimage = single(imread([maskdir,'nucedge_1.tif']));    %open sample to get image size
[height,width] = size(tempimage);
midframe=round(totalframes/2); midframe=midframe+ismember(midframe,frameidx);
sizediffiqr=iqr(wellsss{midframe}(:,4));
x=zeros(1,EF); y=zeros(1,EF);
xdiff=zeros(10,1); ydiff=zeros(10,1);
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
    if f>1
        for i=1:cellnum
            neighbors=find(abs(wellsss{gprev}(:,1)-xx(i))<winrad & abs(wellsss{gprev}(:,2)-yy(i))<winrad);
            sizediff=abs(wellsss{gprev}(neighbors,4)-cellsize(i))/sizediffiqr;
            if numel(neighbors)==1
                sizediff1(i)=sizediff;
                match(i)=sizediff<0.2;
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
save([datadir,'jitter_', shot, '_test.mat'],'x','y');  %later just add this to wellsss
toc(timejittercalc);
cd([codepath,'Processing\']); %return to this directory