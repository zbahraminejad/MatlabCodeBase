function [tracked,newdaughters,nuc_label]=adaptivetrack(xprev,yprev,massprev,initdata,nuc_raw,nuc_label,nucr)
nuc_mask=bwmorph(nuc_label,'remove');
nuc_label_mask=nuc_label.*nuc_mask;
bordermask=zeros(size(nuc_mask));
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=75;
xcur=initdata(:,1);
ycur=initdata(:,2);
masscur=initdata(:,3);
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
nuc_label(mergemask)=newlabels; %label new cells with new IDs
xcur=[xcur;new_center(:,1)]; xcur(mergeid)=[];
ycur=[ycur;new_center(:,2)]; ycur(mergeid)=[];
masscur=[masscur;new_mass];  masscur(mergeid)=[];
numcur=numel(xcur);

%%% match each cell from previous frame to a cell in the current frame %%%
previd=find(~isnan(xprev));
numprev=numel(previd);
newdaughters=ones(numcur,1)*NaN;
xprev=xprev(previd); yprev=yprev(previd);
prevmatch=ones(numprev,2)*NaN;
scores=ones(numprev,2)*NaN;
for i=1:numprev
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    if isempty(neighbors)
        continue;
    end
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    massdiff=(masscur(neighbors)-massprev(i))/massprev(i);
    match=neighbors(massdiff<0.2 && massdiff>-0.2 && dist<winrad);
    if numel(match)==1
        scores(i,1)=dist(match)/winrad+5*massdiff(match);
        prevmatch(i,1)=match;
    elseif numel(match)>1
        %dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
        [tempscore,rank]=sort(dist(match)/winrad+5*massdiff(match));
        scores(i,:)=tempscore([1 2]);
        prevmatch(i,:)=match(rank([1 2]));
    else %check for daughters
        match=neighbors(massdiff>-0.55 && massdiff<-0.40 && dist<20);
        if numel(match)==2
            newdaughters(match(1))=previd(i);
            newdaughters(match(2))=previd(i);
            prevmatch(i,1)=match(1);
            prevmatch=[prevmatch;match(2) NaN];
        end
    end
end

%%% account for every current cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempsort=sort(prevmatch(:,1));
curconflicts=unique(tempsort([diff(tempsort)==0;0]));
for i=curconflicts
    prevcells=find(prevmatch(:,1)==i);
    if numel(prevcells)>2
        prevmatch(prevcells,:)=[NaN NaN];
        continue;
    end
    score1=scores(prevcells(1),:);
    score2=scores(prevcells(2),:);
    config1=nansum([score1(1) score2(2)]);
    config2=nansum([score1(2) score2(1)]);
    if config1<config2
        temp=prevmatch(prevcells(2),:);
        prevmatch(prevcells(2),:)=[temp(2) temp(1)];
    elseif config2>config1
        temp=prevmatch(prevcells(1),:);
        prevmatch(prevcells(1),:)=[temp(2) temp(1)];
    else
        prevmatch(prevcells,:)=ones(2)*NaN;
    end
end