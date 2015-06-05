function bordermask=segmentnucleiandring(set,rlidx,nucmask,bordermask,nucr)
bridgeimage=0;
numpoints=size(set,1);
orderedset=zeros(numpoints,2);
ringid=zeros(numpoints,2);
ringsegid=zeros(numpoints,2);
orderedset(1,:) = set(1,:);

%%% determine clockwise direction from r(1),c(1) %%%%%%%%%%%%%%%%%%
compass=[1,0;1,-1;0,-1;-1,-1;-1,0;-1,1;0,1;1,1];
buffer=zeros(size(set));
x=set(1,1); y=set(1,2);
buffer(:,1)=set(:,1)-x;
buffer(:,2)=set(:,2)-y;
buffer(1,:)=[10000,10000];
adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
if numel(adjidx)~=2
    fprintf('Other than 2 neighbors detected!\n');
end
poles=find(ismember(compass,buffer(adjidx,:),'rows'));
side=compass(round(mean(poles)),:);
nucpres = nucmask(y+side(2),x+side(1));
if nucpres
    pole = poles(1);
else
    pole = poles(2);
end
adjidx=find(buffer(:,1)==compass(pole,1) & buffer(:,2)==compass(pole,2));

%%% order ring coordinates contiguously %%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=2:numpoints-1
    orderedset(pp,:)=set(adjidx,:);
    x=buffer(adjidx,1); y=buffer(adjidx,2);
    if x
                    
        
    buffer(:,1)=buffer(:,1)-x;
    buffer(:,2)=buffer(:,2)-y;
    buffer(adjidx,:)= [10000,10000];
    adjidx=find(abs(buffer(:,1))<=1 & abs(buffer(:,2))<=1);
    if numel(adjidx)~=1
        pole=find(compass(:,1)==x & compass(:,2)==y);
        antipole=pole+4;
        priority=[antipole+1:antipole+7]';
        priority=mod(priority,8); priority(priority==0)=8;
        priority=compass(priority,:);
        poles=find(ismember(priority,buffer(adjidx,:),'rows'));
        x=priority(poles(1),1); y=priority(poles(1),2);
        adjidx=find(buffer(:,1)==x & buffer(:,2)==y);
    end
end
orderedset(numpoints,:)=set(adjidx,:);

%%% calculate tangent angle & detect vertices %%%%%%%%%%%%%%%%%%%%%
offsetshort=round(nucr/4);
offsetlong=2*offsetshort;
gradientoffsetshort=offsetshort; %was offset*2
gradientoffsetlong=offsetlong;

orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
shortdiff=orderedsetoffsetshort-orderedset;
longdiff=orderedsetoffsetlong-orderedset;
shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians
longgrad=atan2(longdiff(:,2),longdiff(:,1));
shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
shortgraddiff=shortgradoffset-shortgrad;
longgraddiff=longgradoffset-longgrad;
shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);
longgraddiff=longgraddiff+2*pi*(longgraddiff<0);
shortgradthresh = pi/6; %was pi/4
longgradthresh = pi/6;
vIdxmasklong=longgraddiff>longgradthresh & longgraddiff<pi;
vIdxmasklong=[zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
vIdxmasklong=imdilate(vIdxmasklong,strel('square',1+nucr));
shortgraddiff(shortgraddiff>=pi)=0;
%vIdxmaskshort=shortgraddiff>shortgradthresh & shortgraddiff<pi;
vIdxmaskshort=shortgraddiff>shortgradthresh;
%vIdxmaskshort=getmax(vIdxmaskshort,shortgraddiff);

vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));
vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
maxmask=zeros(size(vIdxmaskshort));
%vIdx=zeros(length(vIdxobs),1);
for rpc=1:length(vIdxobs)
    %vIdx(rpc)=floor(vIdxobs(rpc).Centroid(2));
    pix=vIdxobs(rpc).PixelIdxList;
    [~,index]=max(shortgraddiff(pix));
    maxmask(pix(index)+offsetshort)=1;
end

vIdxmask=vIdxmasklong & maxmask;
vIdx=find(vIdxmask);

%{
%%% detect missing vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpos=orderedset(vIdx,:);
vIdxmissed=zeros(length(vIdx),1);
tempidx=[1:numpoints]';
pdthresh=nucr*3;
for vc=1:size(vpos,1)
    vdiff=zeros(numpoints,2);
    vdiff(:,1)=orderedset(:,1)-vpos(vc,1);
    vdiff(:,2)=orderedset(:,2)-vpos(vc,2);
    vdist=sqrt(sum(vdiff.^2,2));
    vpd=abs(tempidx-vIdx(vc));
    vpdf=abs(tempidx-(vIdx(vc)+numpoints));
    vpdr=abs(tempidx-(vIdx(vc)-numpoints));
    vIdxreq=vdist<=4 & vpd>pdthresh & vpdf>pdthresh & vpdr>pdthresh;
    vIdxreq=imdilate(vIdxreq,strel('square',nucr+1));
    vIdxreq=find(vIdxreq);
    if ~isempty(vIdxreq) & ~ismember(vIdxreq,vIdx)
        [~,distidx]=sort(vdist(vIdxreq));
        vIdxmissed(vc)=vIdxreq(distidx(1));  %should be distidx(1), but catching errors
    end
end
vIdxmissed=vIdxmissed(vIdxmissed>0);
vIdx=[vIdx;vIdxmissed];
%}

%%% mark candidate vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bridgeimage
    vpos=orderedset(vIdx,:);
    for vc=1:size(vpos,1)
        bordermask(vpos(vc,2),vpos(vc,1))=1;
    end
end
%%% pair and connect vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(vIdx)>=2
    periIdx=vIdx;
    perisize=length(vIdxmask);
    periIdxadj=[periIdx(2:end);perisize+periIdx(1)];
    pairperi=periIdxadj-periIdx;
end
while length(vIdx)>=2
    vpos=orderedset(vIdx,:);
    vposadj=[vpos(2:end,:);vpos(1,:)];
    pair=vposadj-vpos;
    pairdist=sqrt(sum(pair.^2,2));

    [~,ordx]=sort(pairperi./pairdist);
    idx=ordx(end);

    idxadj=idx+1;
    if idxadj>length(vIdx)
        idxadj=1;
    end

    if length(vIdx)<=3
        if length(vIdx)==3
            idxadjadj=idxadj+1;
            if idxadjadj>length(vIdx)
                idxadjadj=1;
            end
            pairperi(idxadj)=pairperi(idxadj)+pairperi(idxadjadj);
            pairperi(idxadjadj)=nucr*pi;
        end
        vlo=pairperi<nucr*pi; %less than half circumference of avg size nucleus
        if sum(vlo)         %either side too short
            break
        end
    end

    [bcx,bcy]=bridge(vpos(idx,:),vpos(idxadj,:));
    for bci=1:length(bcx)
        bordermask(bcy(bci),bcx(bci))=1;
    end

    previdx=idx-1;
    if previdx==0
        previdx=length(vIdx);
    end
    pairperi(previdx)=pairperi(previdx)+length(bcx)-1+pairperi(idxadj);

    vIdx([idx,idxadj])=[];
    pairperi([idx,idxadj])=[];
end
end