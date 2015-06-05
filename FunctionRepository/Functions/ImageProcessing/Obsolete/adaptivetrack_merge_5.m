function [tracedata,curdatatracked,tracking,nuc_label]=adaptivetrack_merge(lastgoodframe,curframe,tracedata,curdata,tracking,nuc_raw,nuc_label,nucr,jitter)
masschangethreshold=0.20; %MCF10A-10xBin1=0.20  BJ5-10xBin2=0.10
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
absjitx=jitter(1); absjity=jitter(2);
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=4*nucr; %farthest jump (pixels)
maxcellnum=size(tracedata,1);
numprevtotal=find(~isnan(tracedata(:,lastgoodframe,1)),1,'last');
xprev=tracedata(1:numprevtotal,lastgoodframe,1); yprev=tracedata(1:numprevtotal,lastgoodframe,2); massprev=tracedata(1:numprevtotal,lastgoodframe,4);
previd=find(~isnan(xprev));
numprevextant=numel(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcurorg=numel(xcur);
ongoingmerge=~isnan(tracking(:,4)) && isnan(tracking(:,5));

%%% attempt to split current or previous merges %%%%%%%%%%%%%%%%%%%%%%%%%%%
bordermask=zeros(size(nuc_mask));
borderflag=zeros(numcurorg,1);
for i=1:numprevextant
    p=previd(i);
    neighbors=find(abs(xcur-xprev(p))<winrad & abs(ycur-yprev(p))<winrad);
    if numel(neighbors)==0
        continue;
    end
    [~,cidx]=min(sqrt((xcur(neighbors)-xprev(p)).^2+(ycur(neighbors)-yprev(p)).^2));
    match=neighbors(cidx);
    massdiff=(masscur(match)-massprev(p))/massprev(p);
    if massdiff>masschangethreshold || ongoingmerge(p)
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
    massdiff=massdiff(candidateidx);
    splitcheck=ongoingmerge(p) & massdiff(1)<-masschangethreshold & massdiff(2)<-masschangethreshold & abs(sum(massdiff)+1)<masschangethreshold;
    daughtercheck=~ongoingmerge(p) & massdiff>-0.55 & massdiff<-0.45 & dist<winrad;
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
        elseif abs(massdiffpre)>masschangethreshold && abs(massdiffpost)>masschangethreshold
            splitorder=isequal(massdiffpre>0,massdiffpost>0);
            if splitorder==1
                cell1post=candidates(1); cell2post=candidates(2);
            else
                cell1post=candidates(2); cell2post=candidates(1);
            end
            tracedata(cell1,mergedframes,:)=(tracedata(cell1post,curframe,:)-tracedata(cell1,premergeframe,:))*(curframe-mergedframes);
            tracedata(cell2,mergedframes,:)=(tracedata(cell2post,curframe,:)-tracedata(cell2,premergeframe,:))*(curframe-mergedframes);
        end
        prevmatch(cell1)=cell1post;
        prevmatch(cell2)=cell2post;
        curmatch(candidates)=curmatch(candidates)+1;
    elseif sum(daughtercheck)==2
        curtracking(candidates,1)=p; %store mothercell
        curmatch(candidates)=curmatch(candidates)+1;
    elseif abs(massdiff(1))<masschangethreshold
        prevmatch(p)=candidates(1);
        curmatch(candidates(1))=p;
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
for c=curconflict'
    prevmatch(prevmatch==c)=NaN;
    curmatch(c)=0;
    curtracking(c,:)=NaN;
end
%%% stop tracking split or non-tracked merged cells %%%%%%%%%%%%%%%%%%%%%%%
tracking(ongoingmerge && isnan(prevmatch),5)=curframe-1;


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