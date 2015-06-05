function [tracedata,curdatatracked,tracking,nuc_label]=adaptivetrack_merge(lastgoodframe,curframe,tracedata,curdata,tracking,nuc_raw,nuc_label,nucr,jitter)
masschangethreshold=0.20; %MCF10A-10xBin1=0.20  BJ5-10xBin2=0.10
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
absjitx=jitter(1); absjity=jitter(2);
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=3*nucr; %prev:4*nucr
maxcellnum=size(tracedata,1);
numprevtotal=find(~isnan(tracedata(:,lastgoodframe,1)),1,'last');
xprev=tracedata(1:numprevtotal,lastgoodframe,1);
yprev=tracedata(1:numprevtotal,lastgoodframe,2);
areaprev=tracedata(1:numprevtotal,lastgoodframe,3);
massprev=tracedata(1:numprevtotal,lastgoodframe,4);
previd=find(~isnan(xprev));
numprevextant=numel(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
ongoingmerge=~isnan(tracking(:,4)) & isnan(tracking(:,5));
numcur=numel(xcur);

%%% match each cell from previous frame to a cell in the current frame %%%%
prevmatch=ones(maxcellnum,1)*NaN;
curmatch=zeros(numcur,1); %number of times matched
curtracking=ones(numcur,3)*NaN; %[mother, mergingcell1, mergingcell2]
newmerge=ones(maxcellnum,1)*NaN;
for i=1:numprevextant
    p=previd(i);
    neighbors=find(abs(xcur-xprev(p))<winrad & abs(ycur-yprev(p))<winrad);
    %%% in case of zero neighbors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(p)).^2+(ycur(neighbors)-yprev(p)).^2);
    massdiff=(masscur(neighbors)-massprev(p))/massprev(p);
    areadiff=(areacur(neighbors)-areaprev(p))/areaprev(p);
    %%% in case of only one neighbor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(neighbors)==1
        if abs(massdiff)<masschangethreshold
            prevmatch(p)=neighbors;
            curmatch(neighbors)=curmatch(neighbors)+1;
        end
        continue;
    end
    %%% for multiple neighbors, find closest 2 neighbors %%%%%%%%%%%%%%%%%%
    [~,didx]=sort(dist);
    candidateidx=[didx(1);didx(2)];
    candidates=neighbors(candidateidx);
    dist=dist(candidateidx);
    massdiff=massdiff(candidateidx); areadiff=areadiff(candidateidx);
    splitcheck=ongoingmerge(p) & massdiff(1)<-masschangethreshold & massdiff(2)<-masschangethreshold & abs(sum(massdiff)+1)<masschangethreshold;
    daughtercheck=~ongoingmerge(p) & massdiff>-0.55 & massdiff<-0.45 & dist<winrad & areadiff<0; %H2B: -55 to -45. NLS: -70 to -30
    if splitcheck
        %tracking(p,5)=curframe-1; %finish tracking merge
        cell1=tracking(p,2); cell2=tracking(p,3);
        mergedframes=tracking(p,4):(curframe-1);
        premergeframe=find(~isnan(tracedata(cell1,:,1)),1,'last');
        mass1pre=tracedata(cell1,premergeframe,4); mass2pre=tracedata(cell2,premergeframe,4);
        massdiffpre=(mass2pre-mass1pre)/mass1pre;
        massdiffpost=(masscur(candidates(2))-masscur(candidates(1)))/masscur(candidates(1));
        if tracking(cell1,1)==tracking(cell2,1)  %sisters
            tracedata([cell1;cell2],mergedframes,:)=tracedata([p;p],mergedframes,:);
            tracedata([cell1;cell2],mergedframes,[3 4])=tracedata([cell1;cell2],mergedframes,[3 4])/2;
            prevmatch(cell1)=candidates(1);
            prevmatch(cell2)=candidates(2);
            curmatch(candidates)=curmatch(candidates)+1;
        elseif abs(massdiffpre)>masschangethreshold && abs(massdiffpost)>masschangethreshold
            splitorder=isequal(massdiffpre>0,massdiffpost>0);
            if splitorder==1
                cell1post=candidates(1); cell2post=candidates(2);
            else
                cell1post=candidates(2); cell2post=candidates(1);
            end
            %tracedata(cell1,mergedframes,:)=(tracedata(cell1post,curframe,:)-tracedata(cell1,premergeframe,:))*(curframe-mergedframes);
            %tracedata(cell2,mergedframes,:)=(tracedata(cell2post,curframe,:)-tracedata(cell2,premergeframe,:))*(curframe-mergedframes);
            prevmatch(cell1)=cell1post;
            prevmatch(cell2)=cell2post;
            curmatch(candidates)=curmatch(candidates)+1;
        end
    elseif sum(daughtercheck)==2
        curtracking(candidates,1)=p; %store mothercell
        curmatch(candidates)=curmatch(candidates)+1;
    elseif abs(massdiff(1))<masschangethreshold
        prevmatch(p)=candidates(1);
        curmatch(candidates(1))=curmatch(candidates(1))+1;
    elseif ~ongoingmerge(p) && massdiff(1)>masschangethreshold
        newmerge(p)=candidates(1);
    end
end
%%% confirm and track new merges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eachmerge=unique(newmerge);
eachmerge(isnan(eachmerge))=[];
for m=eachmerge'
    mergingcells=find(newmerge==m);
    if numel(mergingcells)==2
        massmerge=masscur(m);
        masspre=sum(massprev(mergingcells));
        massdiff=(massmerge-masspre)/masspre;
        if abs(massdiff)<masschangethreshold
            curtracking(m,[2 3])=[mergingcells(1) mergingcells(2)]; %store mergingcells
            curmatch(m)=curmatch(m)+1;
        end
    end
end
%%% remove conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curconflict=find(curmatch>1);
% for c=curconflict'
%     prevmatch(prevmatch==c)=NaN;
%     curmatch(c)=0;
%     curtracking(c,:)=NaN;
% end
for c=curconflict'
    tempprevs=prevmatch==c;
    tempmother=curtracking(c,1);
    if ~isnan(tempmother) %mitoses should be removed
        tempsisters=curtracking(:,1)==tempmother;
        curtracking(tempsisters,1)=NaN;
        curmatch(tempsisters)=0;
    end
    prevmatch(tempprevs)=NaN;
    curmatch(c)=0;
    curtracking(c,:)=NaN;
end
%%% stop tracking split or non-tracked merged cells %%%%%%%%%%%%%%%%%%%%%%%
tracking(ongoingmerge & isnan(prevmatch),5)=curframe-1;
%%% stop tracking chronic merges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chronicmerge=curframe-tracking(:,4)>5 & isnan(tracking(:,5));
tracking(chronicmerge,5)=curframe;

%%% update tracked info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curdatatracked=ones(maxcellnum,4)*NaN;
relabelidx=ones(numcur,1)*NaN;
%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matchprevidx=find(~isnan(prevmatch));
matchcuridx=prevmatch(matchprevidx);
curdatatracked(matchprevidx,:)=[xcur(matchcuridx) ycur(matchcuridx) areacur(matchcuridx) masscur(matchcuridx)];
relabelidx(matchcuridx)=matchprevidx;
%%% add daughters cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curtotal=numprevtotal;
daughtercx=find(~isnan(curtracking(:,1)));
numdaughter=numel(daughtercx);
if numdaughter>0
    daughterpx=(curtotal+1):(curtotal+numdaughter);
    curdatatracked(daughterpx,:)=[xcur(daughtercx) ycur(daughtercx) areacur(daughtercx) masscur(daughtercx)];
    tracking(daughterpx,1)=curtracking(daughtercx,1);
    relabelidx(daughtercx)=daughterpx;
    curtotal=curtotal+numdaughter;
end
%%% add newly merged cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mergecx=find(~isnan(curtracking(:,2)));
nummerge=numel(mergecx);
if nummerge>0
    mergepx=(curtotal+1):(curtotal+nummerge);
    curdatatracked(mergepx,:)=[xcur(mergecx) ycur(mergecx) areacur(mergecx) masscur(mergecx)];
    tracking(mergepx,[2 3])=curtracking(mergecx,[2 3]);
    tracking(mergepx,4)=curframe;
    relabelidx(mergecx)=mergepx;
    curtotal=curtotal+nummerge;
end
%%% add non-tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nonmatchcx=find(curmatch==0);
numnonmatched=numel(nonmatchcx);
if numnonmatched>0
    nonmatchpx=(curtotal+1):(curtotal+numnonmatched);
    curdatatracked(nonmatchpx,:)=[xcur(nonmatchcx) ycur(nonmatchcx) areacur(nonmatchcx) masscur(nonmatchcx)];
    relabelidx(nonmatchcx)=nonmatchpx;
end
%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end