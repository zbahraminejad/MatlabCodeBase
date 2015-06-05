function [tracked,newdaughters,nuc_label]=adaptivetrack(prevdata,curdata,nuc_raw,nuc_label,nucr,jitter,extractmask,reljitx,reljity)
masschangethreshold=0.10; %MCF10A-10xBin1=0.20  BJ5-10xBin2=0.10
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
absjitx=jitter(1); absjity=jitter(2);
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=4*nucr;          %farthest jump (pixels)
mitosisrad=3*nucr;      %farthest jump to daughter (pixels)
xprev=prevdata(:,1); yprev=prevdata(:,2); massprev=prevdata(:,4);
numprevorg=numel(xprev);     %this is the number of traces from beginnning of movie
previd=find(~isnan(xprev));
numprev=numel(previd);
xprev=xprev(previd); yprev=yprev(previd); massprev=massprev(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcurorg=numel(xcur);
mergeid=zeros(numcurorg,1);
bordermask=zeros(size(nuc_mask));
borderflag=zeros(numcurorg,1);
%%% detect merges in current frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numcurorg
	neighbors=find(abs(xprev-xcur(i))<winrad & abs(yprev-ycur(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    %%% find closest neighbor and check mass difference %%%%%%%%%%%%%%%%%%%
    [~,cidx]=min(sqrt((xprev(neighbors)-xcur(i)).^2+(yprev(neighbors)-ycur(i)).^2));
    match=neighbors(cidx);
    massdiff=(masscur(i)-massprev(match))/massprev(match);
    if massdiff>0.1 %possible merge detected
        mergeid(i)=1;
        %%% attempt to segment deflections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%% assign un-merged cells new IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
splitid=find(mergeid & borderflag);
mergemaskorg=ismember(nuc_label,splitid);
mergemask=mergemaskorg & ~bordermask;
mergemask=~bwmorph(~mergemask,'diag');
[newlabels,numnew]=bwlabel(mergemask);
if numnew
    %%% extract features of new cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    new_info=struct2cell(regionprops(newlabels,nuc_raw,'Area','Centroid','MeanIntensity')');
    new_area=squeeze(cell2mat(new_info(1,1,:)));
    new_center=squeeze(cell2mat(new_info(2,1,:)))';
    if numnew==1
        tempx=new_center(1,1)+absjitx;
        tempy=new_center(2,1)+absjity;
        new_center=[tempx tempy];
    else
        new_center(:,1)=new_center(:,1)+absjitx;
        new_center(:,2)=new_center(:,2)+absjity;
    end
    new_density=squeeze(cell2mat(new_info(3,1,:)));
    new_bg=getbackground_H2B(nuc_label,new_center,nuc_raw,nucr,0.25);
    new_density=new_density-new_bg;
    new_mass=new_density.*new_area;
    %%% update with un-merged data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xcur=[xcur;new_center(:,1)]; xcur(splitid)=NaN;
    ycur=[ycur;new_center(:,2)]; ycur(splitid)=NaN;
    areacur=[areacur;new_area];  areacur(splitid)=NaN;
    masscur=[masscur;new_mass];  masscur(splitid)=NaN;
end
numcur=numel(xcur);
%%% update nuc_label with split cell IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label(mergemaskorg>0)=0;
newlabels=newlabels+numcurorg;       %offset so that all label IDs are new
newlabels(newlabels==numcurorg)=0;
nuc_label=nuc_label+newlabels;
%%% match each cell from previous frame to a cell in the current frame %%%%
newdaughters=ones(numcur,1)*NaN;
prevmatch=ones(numprev,1)*NaN;
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    %%% in case of zero neighbors: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    %%% in case of only one neighbor: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(neighbors)==1
        if abs(massdiff)<0.20
            prevmatch(i)=neighbors;
        end
        continue;
    end
    %%% for multiple neighors, find closest 2 neighbors %%%%%%%%%%%%%%%%%%%
    [~,didx]=sort(dist);
    bestidx=didx(1);
    nextidx=didx(2);
    candidateidx=[bestidx;nextidx];
    candidates=neighbors(candidateidx);
    dist=dist(candidateidx);
    massdiff=massdiff(candidateidx);
    %%%%%% check for daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    daughtercheck=massdiff>-0.55 & massdiff<-0.45 & dist<mitosisrad;
    if sum(daughtercheck)==2
        newdaughters(candidates)=previd(i);
        prevmatch=[prevmatch;candidates]; %only used to resolve conflicts
        continue;
    end
    %%%%%% consider closest candidate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if abs(massdiff(1))<masschangethreshold
        prevmatch(i)=candidates(1);
    end
end
%%% resolve conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempsort=sort(prevmatch);
curconflicts=unique(tempsort(logical([diff(tempsort)==0;0])));
for i=curconflicts'
    prevcells=find(prevmatch==i);
    conflictdaughter=prevcells(prevcells>numprev);
    if numel(conflictdaughter)==1 %single daughter should win
        prevmatch(prevcells(prevcells~=conflictdaughter))=NaN;
    elseif numel(conflictdaughter)>1 %must remove the newdaughters assignment
        newdaughters(i)=NaN;
        prevmatch(prevcells)=NaN;
    else
        prevmatch(prevcells)=NaN;
    end
end
%%% update tracked info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prevmatch=prevmatch(1:numprev); %remove daughters
tracked=ones(numprevorg,4)*NaN;
%%%%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempidx=find(~isnan(prevmatch));
matchidxprev=previd(tempidx);
matchidxcur=prevmatch(tempidx);
tracked(matchidxprev,:)=[xcur(matchidxcur) ycur(matchidxcur) areacur(matchidxcur) masscur(matchidxcur)];
%%%%%% add non-tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nontrackedidx=find(~ismember((1:numcur)',prevmatch));
nontracked=[xcur(nontrackedidx) ycur(nontrackedidx) areacur(nontrackedidx) masscur(nontrackedidx)];
mothers=newdaughters(~isnan(newdaughters));
goodnewdaughters=find(~isnan(newdaughters));
newdaughteridx=find(ismember(nontrackedidx,goodnewdaughters));
tracked=[tracked;nontracked];
newdaughters=ones(numel(nontrackedidx),1)*NaN;
newdaughters(newdaughteridx)=mothers;
%%%%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relabelidx=ones(numcur,1)*NaN;
relabelidx(matchidxcur)=matchidxprev;
relabelidx(nontrackedidx)=numprevorg+1:numprevorg+numel(nontrackedidx);
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%NOTE: must include extractmask, jitx, and jity in the arguments
%%%%%% view current cell and its prior neighbors %%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_mask);
dxminprev=max([round(xprev(match)-winrad) 1]); dxmaxprev=min([round(xprev(match)+winrad) width]);
dyminprev=max([round(yprev(match)-winrad) 1]); dymaxprev=min([round(yprev(match)+winrad) height]);
dxmincur=round(dxminprev-reljitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-reljitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-reljity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-reljity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_label,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=extractmask(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=mat2gray(dbmaskprev);
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(dbimage);
%%%%%% view current cell w/ bridge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_mask);
dxminprev=max([round(xprev(match)-winrad) 1]); dxmaxprev=min([round(xprev(match)+winrad) width]);
dyminprev=max([round(yprev(match)-winrad) 1]); dymaxprev=min([round(yprev(match)+winrad) height]);
dxmincur=round(dxminprev-reljitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-reljitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-reljity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-reljity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_label,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbbridgecur=bordermask(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbbridgecur=padarray(dbbridgecur,[dymindiff dxmindiff],'pre');
dbbridgecur=padarray(dbbridgecur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=extractmask(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=dbbridgecur;
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(imresize(dbimage,5));
%%%%%% view previous cell and its future neighbors %%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(nuc_mask);
dxminprev=max([round(xprev(i)-winrad) 1]); dxmaxprev=min([round(xprev(i)+winrad) width]);
dyminprev=max([round(yprev(i)-winrad) 1]); dymaxprev=min([round(yprev(i)+winrad) height]);
dxmincur=round(dxminprev-reljitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-reljitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-reljity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-reljity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_label,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=extractmask(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=mat2gray(dbmaskprev);
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(imresize(dbimage,3));
%}