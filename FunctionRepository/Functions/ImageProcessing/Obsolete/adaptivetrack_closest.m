function [tracked,newdaughters,nuc_label]=adaptivetrack_debug(prevdata,curdata,nuc_raw,nuc_label,nucr,extractmask,jitx,jity)
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
bordermask=zeros(size(nuc_mask));
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=75;     %farthest jump (pixels)
mitosisrad=20; %farthest jump to daughter (pixels)
distcompthresh=nucr; %shortest comparative distance to assume match
xprev=prevdata(:,1); yprev=prevdata(:,2); areaprev=prevdata(:,3); massprev=prevdata(:,4);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcur=numel(xcur);
mergeid=zeros(numcur,1);
borderflag=zeros(numcur,1);
%%% detect merges in current frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numcur
	neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    [~,cidx]=min(sqrt((xprev(neighbors)-xcur(i)).^2+(yprev(neighbors)-ycur(i)).^2));
    match=neighbors(cidx);
    massdiff=(masscur(i)-massprev(match))/massprev(match);
    if massdiff>0.2 %possible merge detected
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
        [bordermask,borderflag(i)]=splitdeflections(orderedset,bordermask,nucr);
    end
end
%%% assign un-merged cells new IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
splitid=find(mergeid & borderflag);
mergemaskorg=ismember(nuc_label,splitid);
mergemask=mergemaskorg & ~bordermask;
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
nuc_label(mergemaskorg>0)=0;
nuc_label=nuc_label+newlabels;
xcur=[xcur;new_center(:,1)]; xcur(splitid)=[];
ycur=[ycur;new_center(:,2)]; ycur(splitid)=[];
masscur=[masscur;new_mass];  masscur(splitid)=[];
numcur=numel(xcur);
%%% match each cell from previous frame to a cell in the current frame %%%
numprevorg=numel(xprev);
previd=find(~isnan(xprev));
numprev=numel(previd);
newdaughters=ones(numcur,1)*NaN;
xprev=xprev(previd); yprev=yprev(previd); massprev=massprev(previd);
prevmatch=ones(numprev,1)*NaN;
distances=zeros(numprev,1);
masschanges=zeros(numprev,1);
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    if numel(neighbors)==1
        if abs(massdiff)<0.2
            prevmatch(i,1)=neighbors;
            distances(i,1)=dist;
            masschanges(i,1)=massdiff;
        end
        continue;
    end
    %%% find closest neighbor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,didx]=sort(dist);
    bestidx=didx(1);
    nextidx=didx(2);
    distcomp=dist(nextidx)-dist(bestidx);
    if distcomp>distcompthresh
        if abs(massdiff(bestidx))<0.2
            prevmatch(i,1)=neighbors(bestidx);
            distances(i,1)=dist(bestidx);
            masschanges(i,1)=massdiff(bestidx);
        end
        continue;
    end
    %%% for close distance differences... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    candidateidx=[bestidx;nextidx];
    candidates=neighbors(candidateidx);
    dist=dist(candidateidx);
    massdiff=massdiff(candidateidx);
    %%%%%% check for daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    daughtercheck=-0.55<massdiff<-0.40 & dist<mitosisrad;
    if sum(daughtercheck)==2
        newdaughters(candidates)=previd(i);
        prevmatch=[prevmatch;candidates]; %only used to resolve conflicts
        distances=[distances;0;0];        %daughter call always wins
        masschanges=[masschanges;0;0];
        continue;
    end
    %%%%%% disqualify inconsistent candidates %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    consistencycheck=find(abs(massdiff)<0.2);
    if numel(consistencycheck)==1
        prevmatch(i,1)=candidates(consistencycheck);
        distances(i,1)=dist(consistencycheck);
        masschanges(i,1)=massdiff(consistencycheck);
    elseif sum(consistencheck)==2
        prevmatch(i,:)=candidates;
        distances(i,:)=dist;
        masschanges(i,:)=massdiff;
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