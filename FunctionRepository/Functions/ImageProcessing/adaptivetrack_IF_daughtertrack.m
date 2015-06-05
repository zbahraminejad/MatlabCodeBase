function [tracked,nuc_label]=adaptivetrack_IF(prevdata,curdata,nuc_label,nucr)
%NOTE: same as adaptivetrack, but mass requirements loosened
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
areaonly=0; %1: Hoechst staining corrupts H2B for mass tracking.
winrad=3*nucr;          %farthest jump (pixels)
masschangethreshold=0.20;
xprev=prevdata(:,1); yprev=prevdata(:,2); areaprev=prevdata(:,3); massprev=prevdata(:,4);
numprevorg=numel(xprev);     %this is the number of traces from beginnning of movie
previd=find(~isnan(xprev));
numprev=numel(previd);
xprev=xprev(previd); yprev=yprev(previd); areaprev=areaprev(previd); massprev=massprev(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcur=numel(xcur);
%%% match each cell from previous frame to a cell in the current frame %%%%
prevmatch=ones(numprev,1)*NaN;
curmatch=zeros(numcur,1); %number of times matched
curtracking=ones(numcur,1)*NaN;
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    %%% in case of zero neighbors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    areadiff=(areacur(neighbors)-areaprev(i))/areaprev(i);
    %%% in case of only one neighbor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if areaonly
        didx=find(dist==min(dist));
        if dist(didx)<nucr && abs(areadiff(didx))<0.1
            matchid=neighbors(didx);
            prevmatch(i)=matchid;
            curmatch(matchid)=curmatch(matchid)+1;
        end
        continue;
    elseif numel(neighbors)==1
        if abs(massdiff)<masschangethreshold
            prevmatch(i)=neighbors;
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
    daughtercheck=massdiff>-0.55 & massdiff<-0.45 & dist<winrad & areadiff<0; %H2B: -55 to -45. NLS: -70 to -30
    if sum(daughtercheck)==2
        curtracking(candidates,1)=i; %store mothercell
        curmatch(candidates)=curmatch(candidates)+1;
    elseif abs(massdiff(1))<masschangethreshold
        prevmatch(i)=candidates(1);
        curmatch(candidates(1))=curmatch(candidates(1))+1;
    end
end
%%% remove conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curconflict=find(curmatch>1);
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
%%%%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracked=ones(numprevorg,4)*NaN;
tempidx=find(~isnan(prevmatch));
matchidxprev=previd(tempidx);
matchidxcur=prevmatch(tempidx);
tracked(matchidxprev,:)=[xcur(matchidxcur) ycur(matchidxcur) areacur(matchidxcur) masscur(matchidxcur)];
%%%%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relabelidx=zeros(numcur,1);
relabelidx(matchidxcur)=matchidxprev;
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end