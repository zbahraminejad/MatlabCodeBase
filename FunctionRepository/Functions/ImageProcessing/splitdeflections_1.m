function [bordermask,bridgeflag]=splitdeflections_1(orderedset,bordermask,nucr)
bridgeflag=0; %returned as 1 if deflections are bridged
nucr=round(nucr/4)*4; %make sure nucr is a multiple of 4
perilength=size(orderedset,1);
%%% detect deflection vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vIdx=getdeflections(orderedset,nucr); %returns boundary indices
%%% count vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vnum=length(vIdx);
if vnum<2
    return; %if less than two vertices are detected, exit function
end
%%% calculate perimeter distance between adjacent vertices %%%%%%%%%%%%%%%%
periIdx=vIdx;
periIdxadj1=[periIdx(2:end);perilength+periIdx(1)];
pairperi1=periIdxadj1-periIdx;
%%% pair and bridge vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while vnum>=2
    skipvertices=0;
    vpos=orderedset(vIdx,:);
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    [bestpair1,ordx1]=sort(pairperi1./pairdist1);
    idx=ordx1(end);
    if idx==vnum
        idxadj=1;
    else
        idxadj=idx+1;
    end
    if vnum>=5
        vposadj2=[vposadj1(2:end,:);vposadj1(1,:)];
        pair2=vposadj2-vpos;
        pairdist2=sqrt(sum(pair2.^2,2));
        pairperi2=[pairperi1(2:end);pairperi1(1)];
        pairperi2=pairperi1+pairperi2;   %perimeter btwn every other vertex
        [bestpair2,ordx2]=sort(pairperi2./pairdist2);
        if bestpair2(end)>bestpair1(end)
            skipvertices=1;
            idx=ordx2(end);
            if idx==vnum
                idxadj=2; idxint=1;
            elseif idx==vnum-1
                idxadj=1; idxint=vnum;
            else
                idxadj=idx+2; idxint=idx+1;
            end
        end
    end

    if vnum<=3
        if vnum==3
            idxadjadj=idxadj+1;
            if idxadjadj>vnum
                idxadjadj=1;
            end
            pairperi1(idxadj)=pairperi1(idxadj)+pairperi1(idxadjadj);
            pairperi1(idxadjadj)=nucr/2*pi;
        end
        vlo=pairperi1<nucr/2*pi;   %less than half circumference of avg size nucleus
        if sum(vlo)             %either side too short
            break
        end
    end

    [bcx,bcy]=bridge(vpos(idx,:),vpos(idxadj,:));
    for bci=1:length(bcx)
        bridgeflag=1;
        bordermask(bcy(bci),bcx(bci))=1;
    end
    
    %%% assign new perimeter distance & remove old vertices %%%%%%%%%%%%%%%
    previdx=idx-1;
    if previdx==0
        previdx=vnum;
    end
    pairperi1(previdx)=pairperi1(previdx)+length(bcx)-1+pairperi1(idxadj);
    if skipvertices==0
        vIdx([idx,idxadj])=[];
        pairperi1([idx,idxadj])=[];
    else
        vIdx([idx,idxint,idxadj])=[];
        pairperi1([idx,idxint,idxadj])=[];
    end
    vnum=length(vIdx);
end

end
%%% debug: visualize deflections on boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
vpos=orderedset(vIdx,:);
for vc=1:size(vpos,1)
    bordermask(vpos(vc,2),vpos(vc,1))=1;
end
%}
