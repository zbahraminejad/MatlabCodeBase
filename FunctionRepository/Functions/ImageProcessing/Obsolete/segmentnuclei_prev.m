function bordermask=segmentnuclei(set,nucmask,bordermask,nucr)
bridgeimage=0;
orderedset=zeros(size(set));
orderedset(1,:) = set(1,:);
numpoints=size(set,1);
%%% determine clockwise direction from r(1),c(1) %%%%%%%%%%%%%%%%%%
x=set(1,1); y=set(1,2);
adj=[x+1,y;x+1,y-1;x,y-1;x-1,y-1;x-1,y;x-1,y+1;x,y+1;x+1,y+1]; %clockwise from pos x axis (wrt y-axis convention)
[Pradj,Ixset]=ismember(adj,set,'rows'); %present in adj, index in set
if sum(Pradj)~=2
    fprintf('Other than 2 neighbors detected!\n');
end
poles=find(Pradj);
side=round(mean(poles));
nucpres = nucmask(adj(side,2),adj(side,1));
if nucpres
    pole = poles(1);
else
    pole = poles(2);
end
current=set(Ixset(pole),:);  %update
set(1,:)=[];                 %remove previous center
%%% order ring coordinates contiguously %%%%%%%%%%%%%%%%%%%%%%%%%%%
for pp=2:numpoints-1
    set((set(:,1)==current(1) & set(:,2)==current(2)),:)=[];    %remove previous center
    orderedset(pp,:)=current;
    x=current(1); y=current(2);
    adj=[x+1,y-1;x,y-1;x-1,y-1;x-1,y;x-1,y+1;x,y+1;x+1,y+1;x+1,y]; %counterclockwise from pos x axis
    [Pradj,Ixset]=ismember(adj,set,'rows'); %present in adj, index in set
    if sum(Pradj)==1
        pole=find(Pradj);
    else
        fprintf('Other than 1 neighbor detected!\n');
        antipole=pole+4;
        priority=[antipole+1:antipole+7]';
        priority=mod(priority,8); priority(priority==0)=8;
        poles=priority(ismember(priority,find(Pradj)));
        pole=poles(1);
    end
    current=set(Ixset(pole),:);  %update
end
orderedset(numpoints,:)=current;
%%% calculate tangent angle & detect vertices %%%%%%%%%%%%%%%%%%%%%
offsetshort=round(nucr/4);
offsetlong=2*offsetshort;
halfoffset=0; %was round(offset/2)
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