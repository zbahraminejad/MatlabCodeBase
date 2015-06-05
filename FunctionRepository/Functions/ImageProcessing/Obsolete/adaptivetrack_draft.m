function [tracked,newdaughters,nuc_label]=adaptivetrack(prevdata,curdata,nuc_raw,nuc_label,nucr)
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
bordermask=zeros(size(nuc_mask));
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=75;
xprev=prevdata(:,1); yprev=prevdata(:,2); massprev=prevdata(:,3);
xcur=curdata(:,1); ycur=curdata(:,2); masscur=curdata(:,3);
numcur=numel(xcur);
mergeid=zeros(numcur,1);
%%% detect merges in current frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numcur
	neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    massdiff=(masscur(i)-massprev(neighbors))./massprev(neighbors);
    match=neighbors(massdiff<0.2);
    if isempty(match) %probable merge detected
        mergeid(i)=1;
        %%% attempt to segment deflections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [r,c]=find(nuc_label_mask==i);
        coorset=[c,r];  %adjust to x-y convention
        [order,status]=orderperimeter([c,r]);
        if status==0    %unable to order perimeter
            fprintf('unable to order perimeter\n');
            continue;
        end
        orderedset=coorset(order,:);
        bordermask=splitdeflections(orderedset,bordermask,nucr);
    end
end
%%% assign un-merged cells new IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mergeid=find(mergeid);
mergemask=ismember(nuc_label,mergeid);
mergemask=mergemask & ~bordermask;
mergemask=~bwmorph(~mergemask,'diag');
newlabels=bwlabel(mergemask);
%%% extract features of new cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_info=struct2cell(regionprops(newlabels,nuc_raw,'Area','Centroid','MeanIntensity')');
new_area=squeeze(cell2mat(new_info(1,1,:)));
new_center=squeeze(cell2mat(new_info(2,1,:)))';
new_density=squeeze(cell2mat(new_info(3,1,:)));
new_bg=getbackground_H2B(nuc_label,new_center,nuc_raw,nucr,0.25);
new_density=new_density-new_bg;
new_mass=new_density.*new_area;
%%% assign new cells new unique IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newlabels=newlabels+numcur;  %offset so that all label IDs are new
newlabels(newlabels==numcur)=0;
%%% update all data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_label(mergemask)=newlabels; %label new cells with new IDs
nuc_label(newlabels>0)=0;
nuc_label=nuc_label+newlabels;
xcur=[xcur;new_center(:,1)]; xcur(mergeid)=[];
ycur=[ycur;new_center(:,2)]; ycur(mergeid)=[];
masscur=[masscur;new_mass];  masscur(mergeid)=[];
numcur=numel(xcur);
%%% match each cell from previous frame to a cell in the current frame %%%
numprevorg=numel(xprev);
previd=find(~isnan(xprev));
numprev=numel(previd);
newdaughters=ones(numcur,1)*NaN;
xprev=xprev(previd); yprev=yprev(previd); massprev=massprev(previd);
prevmatch=ones(numprev,1)*NaN;
scores=zeros(numprev,1);
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    nidx=find(massdiff<0.2 & massdiff>-0.2 & dist<winrad);
    match=neighbors(nidx);
    if numel(match)>=1
        tempscore=dist(nidx)/winrad+5*massdiff(nidx);
        [scores(i),tidx]=min(tempscore);
        prevmatch(i)=match(tidx);
    else %check for daughters
        match=neighbors(massdiff>-0.55 & massdiff<-0.40 & dist<20);
        if numel(match)==2
            newdaughters(match(1))=previd(i);
            newdaughters(match(2))=previd(i);
            prevmatch=[prevmatch;match]; %only used to resolve conflicts
            scores=[scores;0];
        end
    end
end
%%% throw out weaker scores in conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempsort=sort(prevmatch);
curconflicts=unique(tempsort([diff(tempsort)==0;0]));
for i=curconflicts
    prevcells=find(prevmatch==i);
    [~,idx]=min(scores(prevcells));
    prevcells(idx)=[];
    prevmatch(prevcells)=NaN;
end
%%% update tracked info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prevmatch=prevmatch(1:numprev); %remove daughters
tracked=ones(numprevorg,3)*NaN;
relabelidx=ones(numcur,1)*NaN;
%%%%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempidx=find(~isnan(prevmatch));
matchidxorg=previd(tempidx);
matchidx=prevmatch(tempidx);
tracked(matchidxorg,:)=[xcur(matchidx) ycur(matchidx) masscur(matchidx)];
relabelidx(matchidx)=matchidxorg;
%%%%%% add non-tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nontrackedidx=find(~ismember(1:numcur,prevmatch));
nontracked=[xcur(nontrackedidx) ycur(nontrackedidx) masscur(nontrackedidx)];
[mothers,daughteridx]=find(~isnan(newdaughters));
newdaughteridx=ismember(nontrackedidx,daughteridx);
tracked=[tracked;nontracked];
newdaughters=ones(numel(nontrackedidx),1)*NaN;
newdaughters(newdaughteridx)=mothers;
relabelidx(nontrackedidx)=numprevorg+1:numprevorg+1+numel(nontrackedidx);
%%%%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numcur
    nuc_label(nuc_label==i)=relabelidx(i);
end