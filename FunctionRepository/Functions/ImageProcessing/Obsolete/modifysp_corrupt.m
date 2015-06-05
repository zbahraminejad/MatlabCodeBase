function [newsp,trk_rc,corruptlist,leaveout]=modifysp_corrupt(wellsp)
debugmode = 1;
%% create track logic matrix
totalcells = size(wellsp{end},1);
totalframes = size(wellsp,3);
trk_or=zeros(totalcells,totalframes);                   %tracks
siz_or=zeros(size(trk_or));                             %sizes
sig_or=zeros(size(trk_or));                             %signals
for f=1:size(wellsp,3)
    trk_or(1:size(wellsp{f},1),f)=wellsp{f}(:,1)~=0;
    siz_or(1:size(wellsp{f},1),f)=wellsp{f}(:,4);       %nuclear area
    sig_or(1:size(wellsp{f},1),f)=wellsp{f}(:,3);       %median nuclear H2B level
end

%% create track information matrix
trk_df=diff(trk_or,1,2);siz_df=diff(siz_or,1,2);sig_df=diff(sig_or,1,2); % difference between each frame
% set up trk_rc
trk_rc=zeros(size(trk_or,1),5);
%1~5: start frame,split from,end frame,merge to,cell#
% set up siz_vr and sig_vr
siz_vr0=sum(siz_df,2)./sum(siz_or(:,1:end-1),2);siz_vr=prctile(siz_vr0(isfinite(siz_vr0)),[2.5,97.5]);
sig_vr0=sum(sig_df,2)./sum(sig_or(:,1:end-1),2);sig_vr=prctile(sig_vr0(isfinite(sig_vr0)),[2.5,97.5]);

%%%% find each cell's start frame
[r,c]=find(trk_df==1);c=c+1;                            %r=cellid, c=first frame
trk_rc(:,1)=1;trk_rc(r,1)=c;
%%%% find each cell's end frame
[r,c]=find(trk_df==-1);
trk_rc(:,3)=size(trk_or,2);trk_rc(r,3)=c;

%% Assign split/merge to starting/ending cells
id0=1:size(trk_or,1);trk_rc(:,2)=id0;trk_rc(:,4)=id0;trk_rc(:,5)=id0;   %columns 2,4,5 given cell id (i.e. wellsp row)
for f=1:size(wellsp,3)-1
    %% find new cells (new frame) and ending cells (old frame) and minimum distances to nearest cell in other frame
    sid=find(trk_rc(:,1)==f+1);                         %find all tracks that begin at next frame
    smt=ppdist2(wellsp{f+1}(sid,1:2),wellsp{f}(:,1:2)); %dist matrix [next frame new cells, this frame old (all) cells]
    Cs=min(smt,[],2);                                   %find minimum distance for every new cell
    eid=find(trk_rc(:,3)==f);                           %find all tracks that end old frame
    emt=ppdist2(wellsp{f}(eid,1:2),wellsp{f+1}(:,1:2)); %dist matrix [old frame ending cells, new frame all cells
    Ce=min(emt,[],2);                                   %find minimum distance for every ending cell
    f_c=[Cs;Ce];f_th=prctile(f_c,75)+1.5*iqr(f_c);      %calculate outlier level for minimum distances
    %% assign new cell to split from neighbor that: ends, shrinks in size or signal, is closest
    flo=find(wellsp{f}(:,1)>0);fsz=size(flo,1);ssz=size(sid,1);     %flo=cell ids that exist for old frame
    coorfortri=[wellsp{f}(flo,1:2);wellsp{f+1}(sid,1:2)];           %gather coordinates for all cells from old frame and new cells from next frame
    temptri=delaunay(coorfortri(:,1),coorfortri(:,2));
    newsid=fsz+(1:ssz)';                                %returns vector of cellids for new cells which correspond to temptri assignment
    for cc=1:ssz                                        %repeat for each new cell
        [r,dummy]=find(temptri==newsid(cc));            %return rows of temptri that contains new cell's temptri row #
        tempmtx=temptri(r,:);
        tempidx0=unique(tempmtx);tempidx0=tempidx0(tempidx0<=fsz);  %find all unique temptri row# from old frame
        if ~isempty(tempidx0)
            tempidx=flo(tempidx0);                      %get corresponding cell id (i.e. wellsp row #)
            dis=smt(cc,tempidx)';                       %return distances to every unique delaunay-connected cell. *disvec=(dis<f_th);
            dismax=max(dis);
            convec=(wellsp{f+1}(tempidx,1)==0);         %returns 1 for any unique cells whose x coord = 0 in new frame (i.e. ends in old frame)
            %size
            siz_bs=wellsp{f}(tempidx,4);siz_as=siz_bs;          %get nuclear areas of unique cells for old frame
            siz_as(~convec)=wellsp{f+1}(tempidx(~convec),4);    %get nuclear areas of same unique cells for new frame, while convec cells keep prior nuclear area
            sizvec=((siz_as-siz_bs)./siz_bs<=siz_vr(1));        %return 1 for any unique cells whose size increases less than 2.5% percentile calc
            %signal
            sig_bs=wellsp{f}(tempidx,3);sig_as=sig_bs;
            sig_as(~convec)=wellsp{f+1}(tempidx(~convec),3);
            sigvec=((sig_as-sig_bs)./sig_bs<=sig_vr(1));
            %total
            totvec=(convec+(sizvec+sigvec)/2-((dis-dismax)/f_th)*10);maxvec=max(totvec);    %endcell + decreasing area + decreasing signal + closer
            if maxvec>0
                cnd=find(totvec==maxvec);
                [dummy,I]=min(dis(cnd));                %in case two cells have maxvec, find smaller distance
                trk_rc(sid(cc),2)=tempidx(cnd(I));      %assign cell as mother of new cell split
            end
        end
    end
    %% assign ending cell to merge with cell that: begins, increases in size or signal, is closest
    flo=find(wellsp{f+1}(:,1)>0);fsz=size(flo,1);esz=size(eid,1);
    coorfortri=[wellsp{f+1}(flo,1:2);wellsp{f}(eid,1:2)];
    temptri=delaunay(coorfortri(:,1),coorfortri(:,2));
    neweid=fsz+(1:esz)';
    for cc=1:esz
        [r,dummy]=find(temptri==neweid(cc));
        tempmtx=temptri(r,:);
        tempidx0=unique(tempmtx);tempidx0=tempidx0(tempidx0<=fsz);
        if ~isempty(tempidx0)
            tempidx=flo(tempidx0);
            dis=emt(cc,tempidx)';
            dismax=max(dis);
            convec=(tempidx>size(wellsp{f},1));             %any delaunay cells that are new (by way of being added to wellsp)
            %size
            siz_am=wellsp{f+1}(tempidx,4);siz_bm=siz_am;
            siz_bm(~convec)=wellsp{f}(tempidx(~convec),4);
            sizvec=((siz_am-siz_bm)./siz_bm>=siz_vr(2));    %any delaunay cells that increase in size over threshold
            %signal
            sig_am=wellsp{f+1}(tempidx,3);sig_bm=sig_am;
            sig_bm(~convec)=wellsp{f}(tempidx(~convec),3);
            sigvec=((sig_am-sig_bm)./sig_bm>=sig_vr(2));    %same for signal
            %total
            totvec=(convec+(sizvec+sigvec)/2-(dis-dismax)/f_th*10);maxvec=max(totvec);
            if maxvec>0
                cnd=find(totvec==maxvec);
                [dummy,I]=min(dis(cnd));
                trk_rc(eid(cc),4)=tempidx(cnd(I));
            end
        end
    end
end

%% start repairing
orphanlo=zeros(size(trk_rc,1),1);       %any row (cell id) assigned '1' will be removed
mergecorrupted=zeros(size(trk_rc,1),1);         %list of cell id's that merge and will be removed
newsp=wellsp;

%% remove single spots
lo=(trk_rc(:,1)==trk_rc(:,3));
orphanlo(lo)=1;

%% detect missing connection (i.e. cell ends in old frame and new cell starts in next frame, but they are the same cell)
trk_en=trk_rc(:,[3,4,5]);trk_en(:,1)=trk_en(:,1)+1;
trk_st=trk_rc(:,[1,2,5]);
conind=zeros(size(trk_rc,1),2);conlen=0;
%%%% build list of all cell id (ie. trk_rc row #) pairs where an ending
%%%% cell merges with a starting cell or a starting cell splits from an
%%%% ending cell.  Store in matrix conind.
for ff=2:size(trk_or,2)
    lo1=trk_en(:,1)==ff;lo2=trk_st(:,1)==ff;                        %lo1=cells that ended in previous frame. lo2=cells that start in this frame
    trk_ent=trk_en(lo1,:);trk_stt=trk_st(lo2,:);        
    tempmtx=ppdist2(trk_ent(:,2),trk_stt(:,3));                     %diff in cell id # btwn [ending cell merge cell id, starting cell id]
    [conr1,conc1]=find(tempmtx==0);                                 %return indices where diff=0 (i.e. merge cell id = starting cell id)
    tempmtx=ppdist2(trk_ent(:,3),trk_stt(:,2));                     %diff in cell id# btwn [ending cell id, starting cell split cell id]
    [conr2,conc2]=find(tempmtx==0);
    tempind=unique([conr1(:),conc1(:);conr2(:),conc2(:)],'rows');   %unique row-col pairs only. actually, every pair should be duplicated
    tempsz=size(tempind,1);
    if ~isempty(tempind)
        tempind(:,1)=trk_ent(tempind(:,1),3);
        tempind(:,2)=trk_stt(tempind(:,2),3);
        conind(1+conlen:tempsz+conlen,:)=tempind;                   %cell ids of each pair [prev frame, this frame] added to conind
        conlen=conlen+tempsz;
    end
end
if conlen>0
    %%%% Only consider cell pairs as possible missed connections if cell
    %%%% size & signal don't change drastically
    conind=conind(1:conlen,:);
    trk_con=[trk_en(conind(:,1),[1,3]),trk_st(conind(:,2),3)];      %[start frame (actually end+1), end id, start id];
    consiz=[diag(siz_or(trk_con(:,2),trk_con(:,1)-1)),diag(siz_or(trk_con(:,3),trk_con(:,1)))]; %[end cell size, start cell size]
    consig=[diag(sig_or(trk_con(:,2),trk_con(:,1)-1)),diag(sig_or(trk_con(:,3),trk_con(:,1)))]; %[end cell sig, start cell sig]
    consizd=diff(consiz,1,2)./consiz(:,1);                          %size increase as a fraction of last cell size
    consigd=diff(consig,1,2)./consig(:,1);                          %signal increase as a fraction of last cell signal
    consizl=(consizd>=siz_vr(1))&(consizd<=siz_vr(2));              %neither shrinks or grows too much in one step
    consigl=(consigd>=sig_vr(1))&(consigd<=sig_vr(2));              %same for signal
    tlo=consizl&consigl;
    trk_con=trk_con(tlo,:);consigd=consigd(tlo);consizd=consizd(tlo);   %update only for size & sig consistent connections
    
    %%%% check if each connection is exclusive: i.e. no cell merges with
    %%%% two other cells or no cell splits from two cells.
    trk_label=zeros(size(trk_con,1),1);
    condis1=squareform(pdist(trk_con(:,[1,2])));                    %returns distances between each row (pairwise) of trk_con for the 'dimensions' frame and ending cell
    [dis1row,dis1col]=find(condis1==0);dis1lo=(dis1row<dis1col);    %find those where frame & ending cell occurs multiple times
    trk_bad1=[dis1row(dis1lo),dis1col(dis1lo)];                     %only take first instance of pairs (if 1,5 is found, 5,1 will be ignored)
    condis2=squareform(pdist(trk_con(:,[1,3])));                    %repeat for frame & starting cell multiples
    [dis2row,dis2col]=find(condis2==0);dis2lo=(dis2row<dis2col);
    trk_bad2=[dis2row(dis2lo),dis2col(dis2lo)];
    trk_bad=[trk_bad1;trk_bad2];                                    %gather all identical occurances for both ending & starting cells
    if size(trk_bad,1)>0
        %uni_bad=sort(unique(trk_bad(:)));lon_bad=sort(trk_bad(:));
        %if length(uni_bad)<length(lon_bad)
        %   [lon_hist,xi]=hist(lon_bad,uni_bad);worm=xi(lon_hist>1);
        %   for wormc=1:size(worm,1)
        %      [wormr,~]=find(trk_bad==worm(wormc));trk_bad(wormr,:)=[];
        %      trk_label(worm(wormc))=1;
        %   end
        %end
        for cc=1:size(trk_bad,1)
            tempsdf=abs(consigd(trk_bad(cc,:)));                    %2-element column vector of size change differences for each cell
            tempzdf=abs(consizd(trk_bad(cc,:)));
            tempidx=(diff([tempsdf,tempzdf],1,1)>0)-(diff([tempsdf,tempzdf],1,1)<0);    %ex: [1 -1] --> 2nd stronger but smaller
            tix=1+((tempidx(1)*100+tempidx(2))>0);                  %default 1, but if signal change bigger, or same but size change is bigger, then 2
            trk_label(trk_bad(cc,tix))=1;                           %connection with bigger signal change will be removed
        end
    end
    
    %%%% repair missing connection and remove daughter
    trk_con(trk_label==1,:)=[];
    for cc=1:size(trk_con,1)
        x=(trk_rc(:,2)==trk_con(cc,3));trk_rc(x,2)=trk_con(cc,2);   %split from daughter of connection gets reassigned to the connection's orginal cell
        x=(trk_rc(:,4)==trk_con(cc,3));trk_rc(x,4)=trk_con(cc,2);   %merged to daughter reassigned to original
        trk_rc(trk_con(cc,2),[3,4])=trk_rc(trk_con(cc,3),[3,4]);    %original adopts daughter's final frame and merge
        for f=trk_con(cc,1):trk_rc(trk_con(cc,2),3)
            newsp{f}(trk_con(cc,2),:)=newsp{f}(trk_con(cc,3),:);    %newsp (originally wellsss) original data from connection frame to end overwritten with daughter data
        end
        orphanlo(trk_con(cc,3))=1;                                  %add daughter to leave-out list
        trk_con(trk_con(:,2)==trk_con(cc,3),2)=trk_con(cc,2);       %anywhere the daughter is the original in another connection, it gets overwritten as the original of this connection
    end
end

%% correct merge (preparation)
mergelo=(trk_rc(:,4)~=trk_rc(:,5))&(~orphanlo);                     %cells that merge
mermtx=trk_rc(mergelo,3:5);mermtx=sortrows(mermtx,1);               %[last frame, merged cell, cell id] sorted by last frame
mergelo=mermtx(:,3);                                                %merging cell ids sorted by last frame

%% correct merge (start)
for cc=1:size(mermtx,1)
    tempori=trk_rc(mergelo(cc),[3,5]);                              %tempori=[lastframe,cell id] of merging cell
    tempmom=trk_rc(tempori(2),4);                                   %tempmom=[merged cell]
    tempsplit=find((trk_rc(:,2)==tempmom)&(trk_rc(:,1)>tempori(1))&(~orphanlo));    %cells that split from merged cell (tempmom) and are born after merging cell's final frame
    if (trk_rc(tempmom,1)<=tempori(1))                              %merged cell not born after merging cell's final frame
        %possible merge detected
        siz_bm=newsp{tempori(1)}(tempmom,4);                        %merged cell nuclear area before merge
        siz_am=newsp{tempori(1)+1}(tempmom,4);                      %after merge
        sig_bm=newsp{tempori(1)}(tempmom,3);                        %merged cell H2B signal before merge
        sig_am=newsp{tempori(1)+1}(tempmom,3);                      %after merge
        if ((siz_am-siz_bm)/siz_bm>=siz_vr(2))||((sig_am-sig_bm)/sig_bm>=sig_vr(2))     %large increase in size or signal --> true merge
            %true merge detected
            splitlo=0;
            if size(tempsplit,1)>0                                  
                %subsequent split(s) occurs
                tempdau=trk_rc(tempsplit,[1,5]);                    %daughters' [first frame, cell id]
                tempdau=sortrows(tempdau);                          %sort by first frame (then by cell id)
                for ccc=1:size(tempdau,1)
                    if (splitlo==0)
                        tempdaut=tempdau(ccc,:);
                        if (trk_rc(tempmom,3)>=tempdaut(1))         %merged cell's final frame not before daughter's birth
                            siz_bs=newsp{tempdaut(1)-1}(tempmom,4);siz_as=newsp{tempdaut(1)}(tempmom,4);    %merged cell's nuclear size before and after split
                            sig_bs=newsp{tempdaut(1)-1}(tempmom,3);sig_as=newsp{tempdaut(1)}(tempmom,3);    %signal
                            if ((siz_as-siz_bs)/siz_bs<=siz_vr(1))||((sig_as-sig_bs)/sig_bs<=sig_vr(1))     %large decrease in size or signal --> true split
                                % true split
                                splitlo=1;
                                tempdau=tempdaut;
                            end
                        end
                    end
                end
            end
            
            if (splitlo==1)
                % true split after merge
                % fill the gap between ori and dau
                for f=tempori(1)+1:tempdau(1)-1
                    newsp{f}(tempori(2),:)=newsp{f}(tempmom,:);             %merging cell adopts all data from merged cell between merge and split
                end
                % fix tracks after split
                orisig0=newsp{tempori(1)}([tempori(2),tempmom],[3,4]);      %merging cell's last frame [signal,area] for ori & mom
                dausig0=newsp{tempdau(1)}([tempdau(2),tempmom],[3,4]);      %daughter's first frame [signal,area] for dau & mom
                meansig=mean([orisig0;dausig0]);meanmtx=ones(2,1)*meansig;  %identical rows of [avg sig, avg area]
                orisig=orisig0./meanmtx;dausig=dausig0./meanmtx;            %orisig0 & dausig0 normalized over averages
                dismtx=ppdist2(orisig,dausig);                              %distance (sig & area) between each set (ori;mom vs dau;mom)
                
                if (sum(diag(dismtx))<=sum(diag(fliplr(dismtx))))           %ori & dau more similar in size & signal
                    %link dau to ori
                    x=(trk_rc(:,2)==tempdau(2));trk_rc(x,2)=tempori(2);     %future splits off dau redirected to ori
                    x=(trk_rc(:,4)==tempdau(2));trk_rc(x,4)=tempori(2);     %future merges to dau redirected to ori
                    trk_rc(tempori(2),[3,4])=trk_rc(tempdau(2),[3,4]);      %ori final frame & merge to data overwritten w/ dau
                    for f=tempdau(1):trk_rc(tempori(2),3)
                        newsp{f}(tempori(2),:)=newsp{f}(tempdau(2),:);      %adopt all daughter newsp data after split
                    end
                    mergelo(mergelo==tempdau(2))=tempori(2);                %update merging cell list (dau-->ori) for subsequent merge analyses
                else                                                        %ori & mom(after split) more similar in size & signal
                    % link split mom to ori
                    x=(trk_rc(:,1)>=tempdau(1))&(trk_rc(:,2)==tempmom);trk_rc(x,2)=tempori(2);
                    x=(trk_rc(:,3)>=tempdau(1))&(trk_rc(:,4)==tempmom);trk_rc(x,4)=tempori(2);
                    trk_rc(tempori(2),[3,4])=trk_rc(tempmom,[3,4]);
                    for f=tempdau(1):trk_rc(tempori(2),3)
                        newsp{f}(tempori(2),:)=newsp{f}(tempmom,:);
                    end
                    % link dau to merge mom
                    x=(trk_rc(:,2)==tempdau(2));trk_rc(x,2)=tempmom;
                    x=(trk_rc(:,4)==tempdau(2));trk_rc(x,4)=tempmom;
                    trk_rc(tempmom,[3,4])=trk_rc(tempdau(2),[3,4]);
                    for f=tempdau(1):size(trk_or,2);
                        newsp{f}(tempmom,:)=newsp{f}(tempdau(2),:);
                    end
                    mergelo(mergelo==tempmom)=tempori(2);
                    mergelo(mergelo==tempdau(2))=tempmom;
                end
                orphanlo(tempdau(2))=1;                                     %add daughter to leave-out list
            else % no split: add mom to ori
                trk_rc(tempori(2),3)=trk_rc(tempmom,3);                     %overwrite final frame
                if (trk_rc(tempmom,4)==trk_rc(tempmom,5))                   %overwrite merge to
                    trk_rc(tempori(2),4)=trk_rc(tempori(2),5);
                else trk_rc(tempori(2),4)=trk_rc(tempmom,4);
                end
                for f=(tempori(1)+1):trk_rc(tempori(2),3)
                    newsp{f}(tempori(2),:)=newsp{f}(tempmom,:);             %overwrite data after merge
                end
            end
            %%% remove any cells that merge as well as subsequent daughters
            mergecorrupted(tempori(2))=1;
            mergecorrupted(tempmom)=1;
            anylatedaughters = find(trk_rc(:,2)==tempmom & trk_rc(:,1)>tempori(1));
            mergecorrupted(anylatedaughters)=1;
            anylatedaughters = find(trk_rc(:,2)==tempori(2) & trk_rc(:,1)>tempori(1));
            mergecorrupted(anylatedaughters)=1;
        end
    end
end

%% retract 'merge' label for cells that weren't missed connections or true merge-splits or true merges (anyways, validated merges should have no merge label in the end)
lo=(find((trk_rc(:,3)<size(trk_or,2))&(~orphanlo)&(trk_rc(:,4)~=trk_rc(:,5))));
trk_rc(lo,4)=trk_rc(lo,5);

%% remove 'split' label for cells that don't decrease in size or signal dramatically
lo=(find((trk_rc(:,1)>1)&(~orphanlo)&(trk_rc(:,2)~=trk_rc(:,5))));          %find all split events
allmom=trk_rc(lo,[1,2,5]);momsz=size(allmom,1);                             %for each daughter, allmom=[1st frame,mother,cell id]
for cc=1:momsz
    siz_as=newsp{allmom(cc,1)}(allmom(cc,2),4);                             %mother's size after split
    siz_bs=newsp{allmom(cc,1)-1}(allmom(cc,2),4);                           %mother's size before split
    sig_as=newsp{allmom(cc,1)}(allmom(cc,2),3);                             %mother's signal after split
    sig_bs=newsp{allmom(cc,1)-1}(allmom(cc,2),3);                           %mother's signal before split
    if ((siz_as-siz_bs)/siz_bs>siz_vr(1))&&((sig_as-sig_bs)/sig_bs>sig_vr(1))
        % the split is not real; fix daughter's trk_rc
        trk_rc(allmom(cc,3),2)=allmom(cc,3);
    end
end

corrupted = mergecorrupted;
oldlabel = (1:size(orphanlo,1))';
corruptlistorg = oldlabel(corrupted==1);
corruptlist = corruptlistorg;
leaveout = oldlabel(orphanlo==1);
if ~debugmode
    %% remove orphans from newsp
    for f=1:size(newsp,3)
        sz=size(newsp{f},1);
        newsp{f}=newsp{f}(~orphanlo(1:sz),:);               %newsp rows changed
    end
    %%%% start checking from here
    %% relabel trk_rc
    obsoletelabel=oldlabel(orphanlo==1);                    %cell id's (original wellsss row id) of removed cells
    oldlabel=oldlabel(orphanlo==0);                         %cell id's of kept cells
    trk_rc=trk_rc(orphanlo==0,:);                           %trk_rc rows changed
    corruptlist = zeros(size(corrupted,1));
    % correct obsolete labels
    for cc=1:size(obsoletelabel,1)
        obsofind=find(trk_rc(:,2)==obsoletelabel(cc));      %any splits attributed to removed cells are removed (i.e. set to cell id)
        trk_rc(obsofind,2)=obsofind;
        obsofind=find(trk_rc(:,4)==obsoletelabel(cc));      %same for merge
        trk_rc(obsofind,4)=obsofind;
        obsofind=find(trk_rc(:,5)==obsoletelabel(cc));      %same cell id attributed to removed cell renamed new trk_rc row#
        trk_rc(obsofind,5)=obsofind;
    end
    % relabel cell id's with new trk_rc row # (this is consistent with screened newsp
    for cc=1:size(oldlabel,1)
        trk_rc(trk_rc(:,2)==oldlabel(cc),2)=cc;
        trk_rc(trk_rc(:,4)==oldlabel(cc),4)=cc;
        trk_rc(trk_rc(:,5)==oldlabel(cc),5)=cc;
        corruptlist(corruptlistorg==oldlabel(cc))=cc;        %update corruptlist with new row #s (i.e. cell ids)
    end
    leaveout = [];
end
corruptlist = corruptlist(corruptlist>0);