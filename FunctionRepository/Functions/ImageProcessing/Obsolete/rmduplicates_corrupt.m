function [bestsp,best_rc,corruptlist,leaveoutlist]=rmduplicates_corrupt(newsp,trk_rc,corrupted,nucr,leaveout)
debugmode = 1;
%% logic for duplicates
bestsp=newsp;
best_rc=trk_rc;
duplo=zeros(size(best_rc,1),1);
combi=unique(trk_rc(:,[1,3]),'rows');scombi=size(combi,1);                      %also sorted, first by start then by stop
roundnum=(scombi^2+scombi)/2;
counter=0;tempper0=0;

%% main nested loop
for f1=1:scombi
    %% find tracks starting f1 and ending f1
    TRK1lo=find((best_rc(:,1)==combi(f1,1))&(best_rc(:,3)==combi(f1,2))&(~ismember(best_rc(:,5),leaveout)));       %all cells that run from given start and stop
    sz1=size(TRK1lo,1);pt1=ceil(sz1/1000);
    for f2=f1:scombi
        counter=counter+1;
        tempper=round(counter/roundnum*100);
        if (mod(tempper,10)==0)&&(mod(tempper0,10)~=0)
            disp(['Remove duplicates......',num2str(tempper),'%'])
        end
        tempper0=tempper;
        %% find tracks starting f2 and ending f2
        TRK2lo=find((best_rc(:,1)==combi(f2,1))&(best_rc(:,3)==combi(f2,2))&(~ismember(best_rc(:,5),leaveout)));   %all cells that start and stop in later tracks
        sz2=size(TRK2lo,1);pt2=ceil(sz2/1000);
        %% find track pairs
        for row=1:pt1
            y1=1+1000*(row-1);y2=min([1000*row,sz1]);
            for col=1:pt2
                x1=1+1000*(col-1);x2=min([1000*col,sz2]);
                %% get the distance matrix of each frame
                frnum=0;
                for f=combi(f2,1):min(combi([f1,f2],2))                         %later start to earlier stop
                    %%
                    if frnum==0;
                        DIS0=ppdist2(bestsp{f}(TRK1lo(y1:y2),1:2),bestsp{f}(TRK2lo(x1:x2),1:2));
                        [tempr0,tempc0]=find(DIS0<6*nucr);
                        tempr=tempr0+y1-1;tempc=tempc0+x1-1;
                        if (f1==f2)
                            templo=(tempr<tempc);                               %pairs will be duplicated, so consider 1st instanct
                            tempr=tempr(templo);tempc=tempc(templo);
                        end
                        DIS=diag(DIS0(tempr-y1+1,tempc-x1+1));                  %get distance between each cell
                        dis0=(DIS==0);                                          %1 wherever distance is zero
                        r=TRK1lo(tempr);c=TRK2lo(tempc);                        %get original cell ids
                        duplotemp=((duplo(r)==0)&(duplo(c)==0));                %only consider cells that haven't been called duplicates yet
                        DIS=DIS(duplotemp);dis0=dis0(duplotemp);                %DIS=distance between cells; dis0=non duplicate same coordinates
                        r=r(duplotemp);c=c(duplotemp);                          %cells that share location
                    elseif ~isempty(r)
                        disvec=zeros(size(dis0,1),1);                           
                        for cc=1:size(r)
                            disvec(cc)=pdist(bestsp{f}([r(cc),c(cc)],1:2));     %calc distances between previously identified cells
                        end
                        dislo=(disvec<4*nucr);                                  %must stay within 4 nuclear radii of each other
                        disvec=disvec(dislo);DIS=DIS(dislo);dis0=dis0(dislo);   
                        r=r(dislo);c=c(dislo);
                        DIS=DIS+disvec;                                         %cumulative distance between each pair kept
                        dis0=dis0+(disvec==0);                                  %cumulative # of frames cells occupy same point
                    end
                    frnum=frnum+1;
                end
                %% find track pairs
                DIS=DIS/frnum;dis0=dis0/frnum;                                  %DIS=avg dist btwn pairs. dis0=fraction of time they were segemented as the same cell
                DISLO=find(((DIS<1.5*nucr)&(dis0>0.5)));                        %must be close enough & frequently segmented as same cell
                r=r(DISLO);c=c(DISLO);
                
                %% correct track pairs
                if ~isempty(r)
                    %label duplicates
                    duplo(c)=1;
                    %relabel best_rc
                    best_rc(r,3)=max([best_rc(r,3),best_rc(c,3)],[],2);         %whichever cell had latest last frame
                    for cc=1:size(r)
                        best_rc(best_rc(:,2)==c(cc),2)=r(cc);                   %adopt any daughters of duplicate
                    end
                    %average coordinates and values in bestsp
                    for f=combi(f2,1):min(combi([f1,f2],2))
                        mt1=bestsp{f}(r,:);mt2=bestsp{f}(c,:);
                        bestsp{f}(r,:)=(mt1+mt2)/2;
                    end
                    %remainder of track2 copied to track1
                    if combi(f1,2)<combi(f2,2)
                        for f=(combi(f1,2)+1):combi(f2,2)
                            bestsp{f}(r,:)=bestsp{f}(c,:);
                        end
                    end
                end
            end
        end
    end
end

corruptlist = corrupted;
oldlabel=(1:size(duplo,1))';
obsoletelabel=oldlabel(duplo==1);
leaveoutlist = [leaveout;obsoletelabel];
if ~debugmode
    %% remove duplicates and relabel best_rc
    oldlabel=oldlabel(duplo==0);
    best_rc(duplo==1,:)=[];
    corruptlist = zeros(size(corrupted,1));
    % correct obsolete labels
    for cc=1:size(obsoletelabel,1)
        obsofind=find(best_rc(:,2)==obsoletelabel(cc));
        best_rc(obsofind,2)=obsofind;
        obsofind=find(best_rc(:,4)==obsoletelabel(cc));
        best_rc(obsofind,4)=obsofind;
        obsofind=find(best_rc(:,5)==obsoletelabel(cc));
        best_rc(obsofind,5)=obsofind;
    end
    % relabel remains
    for cc=1:size(oldlabel,1)
        best_rc(best_rc(:,2)==oldlabel(cc),2)=cc;
        best_rc(best_rc(:,4)==oldlabel(cc),4)=cc;
        best_rc(best_rc(:,5)==oldlabel(cc),5)=cc;
        corruptlist(corrupted==oldlabel(cc))=cc;        %update corruptlist with new row #s (i.e. cell ids)
    end
    corruptlist = corruptlist(corruptlist>0);
    %% remove duplicates in bestsp
    for f=1:size(bestsp,3)
        fsz=size(bestsp{f},1);
        bestsp{f}(duplo(1:fsz)==1,:)=[];
    end
    leaveoutlist = [];
end