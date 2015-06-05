function reportmistracking(rawdir,nucr,tracedata,genealogy,jitters,tracestats)
rawdummy=single(imread([rawdir,'CFP_1.tif']));
[maxy,maxx]=size(rawdummy);
[numtraces,numframes,~]=size(tracedata);
fin=find(tracestats(:,2)<numframes-1);
nodaughters=zeros(numel(fin),1);
for i=1:numel(fin)
    nodaughters(i)=isempty(find(genealogy==fin(i),1));
end
fin=fin(logical(nodaughters));
finend=tracestats(fin,2);
orphan=find(tracestats(:,1)>1 & isnan(genealogy));
ostart=tracestats(orphan,1);
numorph=numel(orphan);
omass=ones(numorph,1)*NaN; ox=omass; oy=omass;
for i=1:numel(orphan)
    omass(i)=tracedata(orphan(i),ostart(i),4); ox(i)=tracedata(orphan(i),ostart(i),1); oy(i)=tracedata(orphan(i),ostart(i),2);
end
neighborrad=3*nucr;
borderrad=4*nucr;
breakcount=0;
for i=1:numel(fin)
    fmass=tracedata(fin(i),finend(i),4); fx=tracedata(fin(i),finend(i),1); fy=tracedata(fin(i),finend(i),2);
    fprintf('%0.0f: F%0.0f X%0.0f Y%0.0f M%0.0f\n',i,finend(i),fx-jitters(finend(i),1),fy-jitters(finend(i),2),fmass);
    if fx<borderrad || maxx-fx<borderrad || fy<borderrad || maxy-fy<borderrad %ignore cells close to border
        continue;
    end
    breakcount=breakcount+1;
    c=find(ostart-finend(i)==1 & abs(ox-fx)<neighborrad & abs(oy-fy)<neighborrad);
    if isempty(c)
        continue;
    elseif numel(c)==1
        c=[c;c];
    end
    operc1=100*(omass(c(1))-fmass)/fmass; operc2=100*(omass(c(2))-fmass)/fmass;
    ojx=jitters(ostart(c(1)),1); ojy=jitters(ostart(c(1)),2);
    fprintf('\t\tX%0.0f Y%0.0f M%0.0f  X%0.0f Y%0.0f M%0.0f\n',ox(c(1))-ojx,oy(c(1))-ojy,operc1,ox(c(2))-ojx,oy(c(2))-ojy,operc2);
end
fprintf('%0.0f broken traces out of %0.0f total traces\n',breakcount,numtraces);
minframe=min(tracestats(:,1));maxframe=max(tracestats(:,2));
hist(finend,minframe:maxframe);xlim([minframe maxframe]);