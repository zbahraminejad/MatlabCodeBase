function [alldata,curdatafinal,tracking,nuc_label]=adaptivetrack_merge(lastgoodframe,curframe,alldata,curdata,tracking,nuc_raw,nuc_label,nucr,jitter)
masschangethreshold=0.20; %MCF10A-10xBin1=0.20  BJ5-10xBin2=0.10
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
absjitx=jitter(1); absjity=jitter(2);
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=4*nucr;          %farthest jump (pixels)
xprev=alldata{lastgoodframe}(:,1); yprev=alldata{lastgoodframe}(:,2); massprev=alldata{lastgoodframe}(:,4);
numprevorg=numel(xprev);     %this is the number of traces from beginnning of movie
previd=find(~isnan(xprev));
numprev=numel(previd);
xprev=xprev(previd); yprev=yprev(previd); massprev=massprev(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcurorg=numel(xcur);
%%% attempt to split current or previous merges %%%%%%%%%%%%%%%%%%%%%%%%%%%
bordermask=zeros(size(nuc_mask));
borderflag=zeros(numcurorg,1);
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    if numel(neighbors)==0
        continue;
    end
    [~,cidx]=min(sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2));
    match=neighbors(cidx);
    massdiff=(masscur(match)-massprev(i))/massprev(i);
    if massdiff>masschangethreshold || tracking(previd(i),2)==1
        [bordermask,borderflag(match)]=attemptsplit(match,nuc_label_mask,bordermask,nucr);
    end
end
%%% assign new IDs to un-merged cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
splitid=find(borderflag);
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
        new_center=[new_center(1,1) new_center(2,1)];
    end
    new_density=squeeze(cell2mat(new_info(3,1,:)));
    new_bg=getbackground_H2B(nuc_label,new_center,nuc_raw,nucr,0.25);
    new_density=new_density-new_bg;
    new_mass=new_density.*new_area;
    new_center(:,1)=new_center(:,1)+absjitx;
    new_center(:,2)=new_center(:,2)+absjity;
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
%newdaughters=ones(numcur,1)*NaN;

%newsplits=ones(numcur,1)*NaN;
prevmatch=ones(numprev,1)*NaN;
curmatch=zeros(numcur,1);
curtracking=ones(numcur,1)*NaN;
%removemerge=zeros(numprevorg,1);
%newmerge=ones(numprev,1)*NaN;
for i=1:numprev
    p=previd(i);
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    %%% in case of zero neighbors: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        %removemerge(p)=tracking(p,2)==1; %stop tracking merge
        tracking(p,2)=0;
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    %%% in case of only one neighbor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(neighbors)==1
        if abs(massdiff)<masschangethreshold
            prevmatch(i)=neighbors;
            curmatch(neighbors)=curmatch(neighbors)+1;
        else
            %removemerge(p)=tracking(p,2)==1;
            tracking(p,2)=0;
        end
        continue;
    end
    %%% for multiple neighbors, find closest 2 neighbors %%%%%%%%%%%%%%%%%%
    [~,didx]=sort(dist);
    candidateidx=[didx(1);didx(2)];
    candidates=neighbors(candidateidx);
    dist=dist(candidateidx);
    massdiff=massdiff(candidateidx);
    splitcheck=tracking(p,2)==1 & massdiff(1)<-masschangethreshold & massdiff(2)<-masschangethreshold & abs(sum(massdiff)+1)<masschangethreshold;
    daughtercheck=tracking(p,2)==0 & massdiff>-0.55 & massdiff<-0.45 & dist<winrad;
    if splitcheck
        %newsplits(candidates)=p;
        %prevmatch=[prevmatch;candidates];
        cell1=tracking(p,3); cell2=tracking(p,4);
        if tracking(cell1,1)==tracking(cell2,1)  %sisters
            alldata;
        end
        tracking(p,2)=0;
        tracking(p,5)=curframe;
        
        curmatch(candidates)=curmatch(candidates)+1;
    elseif sum(daughtercheck)==2
        newdaughters(candidates)=p;
        prevmatch=[prevmatch;candidates];
    elseif abs(massdiff(1))<masschangethreshold
        prevmatch(i)=candidates(1);
    elseif massdiff(1)>masschangethreshold
        newmerge(i)=candidates(1);
    end
end
%%% resolve conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempsort=sort(prevmatch);
curconflicts=unique(tempsort(logical([diff(tempsort)==0;0])));
for i=curconflicts'
    prevcells=find(prevmatch==i);
    conflictpriority=prevcells(prevcells>numprev);
    if numel(conflictpriority)==1 %single daughter or split wins
        prevmatch(prevcells(prevcells~=conflictpriority))=NaN;
    elseif numel(conflictpriority)>1 %remove split or daughter assignment
        newsplits(i)=NaN;
        newdaughters(i)=NaN;
        prevmatch(prevcells)=NaN;
    else
        prevmatch(prevcells)=NaN;
    end
end
%%% update tracked info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prevmatch=prevmatch(1:numprev); %remove daughters
curdatatracked=ones(numprevorg,4)*NaN;
%%%%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempidx=find(~isnan(prevmatch));
matchidxprev=previd(tempidx);
matchidxcur=prevmatch(tempidx);
curdatatracked(matchidxprev,:)=[xcur(matchidxcur) ycur(matchidxcur) areacur(matchidxcur) masscur(matchidxcur)];
%%%%%% add non-tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nontrackedidx=find(~ismember((1:numcur)',prevmatch));
curdatanottracked=[xcur(nontrackedidx) ycur(nontrackedidx) areacur(nontrackedidx) masscur(nontrackedidx)];
mothers=newdaughters(~isnan(newdaughters));
goodnewdaughters=find(~isnan(newdaughters));
newdaughteridx=find(ismember(nontrackedidx,goodnewdaughters));
curdatafinal=[curdatatracked;curdatanottracked];
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