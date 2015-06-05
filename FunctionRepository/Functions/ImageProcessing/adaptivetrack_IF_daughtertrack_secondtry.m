function [tracked,nuc_label]=adaptivetrack_IF(prevdata,curdata,nuc_label,nucr,debugpackage)
%NOTE: track daughters (by tracking mass)
H2BvsNLS=1;
if H2BvsNLS==1
    daughtermassmin=-0.60; %timelapse:-0.55 IF:-0.50
    daughtermassmax=-0.40; %timelapse:-0.45 IF:-0.40
    daughterareamin=-0.50; %IF:-0.45
    daughterareamax=-0.15; %IF:-0.25
    daughterintmin=-0.30; %IF:-0.30
    daughterintmax=0.10; %IF:0.05
elseif H2BvsNLS==2
    daughtermassmin=-0.70; %timelapse:-0.70
    daughtermassmax=-0.30; %timelapse:-0.30
end
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=3*nucr;          %farthest jump (pixels)
masschangethreshold=0.25;
xprev=prevdata(:,1); yprev=prevdata(:,2); areaprev=prevdata(:,3); massprev=prevdata(:,4); intprev=massprev./areaprev;
numprevorg=numel(xprev);     %this is the number of traces from beginning of movie
previd=find(~isnan(xprev));
numprev=numel(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4); intcur=masscur./areacur;
numcur=numel(xcur);
%%% match each cell from previous frame to a cell in the current frame %%%%
prevmatch=ones(numprevorg,1)*NaN;
curmatch=zeros(numcur,1); %number of times matched
curtracking=ones(numcur,1)*NaN; %mother
prevtemp=ones(numprevorg,2)*NaN; massdiffrec=prevtemp; areadiffrec=prevtemp; intdiffrec=prevtemp;
for i=1:numprev
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
    intdiff=(intcur(neighbors)-intprev(p))/intprev(p);
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
    massdiff=massdiff(candidateidx); areadiff=areadiff(candidateidx); intdiff=intdiff(candidateidx);
    %daughtercheck=massdiff>daughtermin & massdiff<daughtermax & dist<winrad & areadiff<0; %H2B: -0.55 to -0.45; NLS: -0.70 to -0.30
    daughtercheck=massdiff>daughtermassmin & massdiff<daughtermassmax & dist<winrad & areadiff>daughterareamin & areadiff<daughterareamax & intdiff>daughterintmin & intdiff<daughterintmax; %H2B: -0.55 to -0.45; NLS: -0.70 to -0.30
    if sum(daughtercheck)==2
        curtracking(candidates)=p; %store mothercell
        curmatch(candidates)=curmatch(candidates)+1;
    elseif abs(massdiff(1))<masschangethreshold
        prevmatch(p)=candidates(1);
        curmatch(candidates(1))=curmatch(candidates(1))+1;
    end
    massdiffrec(p,:)=massdiff; areadiffrec(p,:)=areadiff; intdiffrec(p,:)=intdiff;
end
%%% remove conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unique(curtracking(~isnan(curtracking)))
%prevjitter=debugpackage{4}; prevjitx=prevjitter(1); prevjity=prevjitter(2);
%tp=3767; massdiffrec(tp,:),areadiffrec(tp,:),intdiffrec(tp,:)
curconflict=find(curmatch>1);
for c=curconflict'
    tempprevs=find(prevmatch==c);
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
%tempidx=find(~isnan(prevmatch));
%matchidxprev=previd(tempidx);
%matchidxcur=prevmatch(tempidx);
matchidxprev=find(~isnan(prevmatch));
matchidxcur=prevmatch(matchidxprev);
tracked(matchidxprev,:)=[xcur(matchidxcur) ycur(matchidxcur) areacur(matchidxcur) masscur(matchidxcur)];
%%%%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relabelidx=zeros(numcur,1);
relabelidx(matchidxcur)=matchidxprev;
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
previdx=3767;
prevcellandcurrentneighborsIF(debugpackage,winrad,xprev,yprev,previdx);
%}