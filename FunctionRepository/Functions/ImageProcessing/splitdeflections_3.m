function [bordermask,bridgeflag]=splitdeflections_3(orderedset,bordermask,nucr)
%orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)]; %only if using bwboundaries
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
    %%% Determine adjacent vertices that define the highest curvature %%%%%
    vposadj1=[vpos(2:end,:);vpos(1,:)];
    pair1=vposadj1-vpos;
    pairdist1=sqrt(sum(pair1.^2,2));
    curvature1=pairperi1./pairdist1;
    [bestcurve1,curve1idx]=sort(curvature1);
    bestcurveidx=curve1idx(end);
    if bestcurveidx==vnum
        bestcurveidxadj=1;
    else
        bestcurveidxadj=bestcurveidx+1;
    end
    %%% If skipping 1 vertex gives a better curve, take that instead %%%%%%
    if vnum>=5
        vposadj2=[vposadj1(2:end,:);vposadj1(1,:)];
        pair2=vposadj2-vpos;
        pairdist2=sqrt(sum(pair2.^2,2));
        pairperi2=[pairperi1(2:end);pairperi1(1)];
        pairperi2=pairperi1+pairperi2;   %perimeter btwn every other vertex
        curvature2=pairperi2./pairdist2;
        [bestcurve2,curve2idx]=sort(curvature2);
        if bestcurve2(end)>bestcurve1(end)
            skipvertices=1;
            bestcurveidx=curve2idx(end);
            if bestcurveidx==vnum
                bestcurveidxadj=2; idxint=1;
            elseif bestcurveidx==vnum-1
                bestcurveidxadj=1; idxint=vnum;
            else
                bestcurveidxadj=bestcurveidx+2; idxint=bestcurveidx+1;
            end
        end
    end
    %%% If only 3 vertices, assume the 2nd vertex following the vertex
    %%% initiating the best curve is wrong and remove it.
    if vnum<=3
        if vnum==3
            bestcurveidxadj2=bestcurveidxadj+1;
            if bestcurveidxadj2>vnum
                bestcurveidxadj2=1;
            end
            % Given vertex 1-2 gives best curvature, assume vertex 3 is
            % wrong and ignore it (i.e. perimeter 2-1 = p(2-3)+p(3-1))
            pairperi1(bestcurveidxadj)=pairperi1(bestcurveidxadj)+pairperi1(bestcurveidxadj2);
            pairperi1(bestcurveidxadj2)=nucr/2*pi;
        end
        % If either perimeter is less than cutoff, stop segmenting
        vlo=pairperi1<nucr/2*pi;
        if sum(vlo)
            break
        end
    end
    %%% If this point is reached, a split will be performed, so mark it. %%
    bridgeflag=1;
    %%% Bridge the vertices defining the best curvature %%%%%%%%%%%%%%%%%%%
    [bx,by]=bridge(vpos(bestcurveidx,:),vpos(bestcurveidxadj,:));
    for bci=1:length(bx)
        %bridgeflag=1;
        bordermask(by(bci),bx(bci))=1;
    end
    %%% assign new perimeter distances & remove old vertices %%%%%%%%%%%%%%
    previdx=bestcurveidx-1;
    if previdx==0
        previdx=vnum;
    end
    % Given vertex 3-4 gave best curvature, and is now bridged, define the
    % perimeter from vertex 2 to 5: p(2-5)=p(2-3)+bridge+p(4-5).
    pairperi1(previdx)=pairperi1(previdx)+length(bx)-1+pairperi1(bestcurveidxadj);
    % Remove the vertices and perimeters of the vertices defining (or
    % intervening) the best curve.
    if skipvertices==0
        vIdx([bestcurveidx,bestcurveidxadj])=[];
        pairperi1([bestcurveidx,bestcurveidxadj])=[];
    else
        vIdx([bestcurveidx,idxint,bestcurveidxadj])=[];
        pairperi1([bestcurveidx,idxint,bestcurveidxadj])=[];
    end
    vnum=length(vIdx);
end
%keyboard;
end
%%% debug: visualize deflections on boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% vpos=orderedset(vIdx,:);
% for vc=1:size(vpos,1)
%     bordermask(vpos(vc,2),vpos(vc,1))=1;
% end

tempbridgeimage=zeros(size(scmask));
for bci=1:length(bcx)
   tempbridgeimage(bcy(bci),bcx(bci))=1;
end
tempimage=mat2gray(scmask);
tempimage(:,:,2)=tempbridgeimage;
tempimage(:,:,3)=0;
imshow(tempimage);
%}
