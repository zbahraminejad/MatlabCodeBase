function [tracked,nuc_label]=adaptivetrack_IF(prevdata,curdata,nuc_label,nucr,areaonly)
%NOTE: same as adaptivetrack, but mass requirements loosened
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    didx=find(dist==min(dist));
    matchid=neighbors(didx);
    if (areaonly && dist(didx)<nucr && abs(areadiff(didx))<0.2) || ...
            (~areaonly && abs(massdiff(didx))<masschangethreshold)
        prevmatch(i)=matchid;
        curmatch(matchid)=curmatch(matchid)+1;
    end
end
%%% remove conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curconflict=find(curmatch>1);
for c=curconflict'
    tempprevs=prevmatch==c;
    prevmatch(tempprevs)=NaN;
    curmatch(c)=0;
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