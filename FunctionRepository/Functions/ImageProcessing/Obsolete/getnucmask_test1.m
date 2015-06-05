function DAs_pad=getnucmask(DAs_bs,nucr)

[DAs_ma,th1,th2,bg]=ThreshImage_MC_test(DAs_bs,0,25);  %10 bins higher for stricter threshold was log:60 or abs:30

DAs_ma=imfill(DAs_ma,'holes');
DAs_pad=padarray(DAs_ma,[1 1]);
DAs_bs=padarray(DAs_bs,[1 1]);
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));

%{
DAs_ma=DAs_pad(2:end-1,2:end-1);
DAs_bs(DAs_ma==0)=0;
cytogradients=edge(DAs_bs,'canny',[0.0001 0.003],1); %[.0063 .0156]
cytogradients=bwmorph(cytogradients,'thin',Inf);
cytogradients=imclose(cytogradients,strel('disk',2,0));
%cytogradients=imfill(cytogradients,'holes');
%cytogradients=imopen(cytogradients,strel('disk',1,0));
%}

%DAs_pad=padarray(DAs_ma,[1 1]);
%DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));  %this would remove stuff smaller than 25% of normal cell size.  Also removes spurs at edges & protruding points (but not intruding points, I think).
fig1=figure; fig2=figure; fig3=figure;
[nuclabels,obnum]=bwlabel(DAs_pad);
for ci=1:obnum       %ci = 1,2,12
    pixints=DAs_bs(find(nuclabels==ci));
    [ksn,ksx]=ksdensity(pixints,'width',0.05);
    
    
    pixidx=find(nuclabels==ci);
    tmap=zeros(size(DAs_bs));
    tmap(pixidx)=mat2gray(DAs_bs(pixidx));
    [tmapr,tmapc]=find(tmap>0);
    tmap=tmap(min(tmapr):max(tmapr),min(tmapc):max(tmapc));
    set(0,'currentfigure',fig1);
    imshow(imresize(tmap,5));
    set(0,'currentfigure',fig2);
    plot(ksx,ksn);
    set(0,'currentfigure',fig3);
    hist(pixints,100);
    %}
    
    %{
    firstpeak=find(diff(ksn)<0,1,'first');
    %if ksx(firstpeak)>th2
    %    continue
    %end
    nextvalley=find(diff(ksn)>0);
    nextvalley(nextvalley<=firstpeak)=[];
    if isempty(nextvalley)
        continue
    end
    nextvalley=nextvalley(1);
    pixidx=find(nuclabels==ci);
    pixtodelete=pixidx(DAs_bs(pixidx)<ksx(nextvalley));
    if numel(pixtodelete)*2<numel(pixidx)
        DAs_pad(pixtodelete)=0;
    end
    %}
end

DAs_pad=imfill(DAs_pad,'holes');
DAs_pad=imopen(DAs_pad,strel('disk',nucr/2,0));
DAs_pad=~bwmorph(~DAs_pad,'diag');         %break connections
DAs_pad=~bwmorph(~DAs_pad,'bridge');      %break connections
nucedges=bwmorph(DAs_pad,'remove');
[ringlabels,obnum]=bwlabel(nucedges);
bordermask=zeros(size(DAs_pad));
%%% border detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci=1:obnum
    [r,c]=find(ringlabels==ci);
    coorset=[c,r];  %adjust to x-y convention
    %{
    if ci==16   %32
        fprintf('ci = %0.0f\n',ci);
        pixints=DAs_bs(find(nuclabels==ci));
        pixidx=find(nuclabels==ci);
        tmap=zeros(size(DAs_pad));
        tmap(pixidx)=1;
        imshow(tmap);
        
        pixselect=pixidx(DAs_bs(pixidx)>th2);
        tmap=zeros(size(DAs_pad));
        tmap(pixselect)=1;
        imshow(tmap);
    end
    %}
    order=orderperimeter_nuc(coorset);
    orderedset=coorset(order,:);
    bordermask=segmentnuclei(orderedset,bordermask,nucr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DAs_pad=DAs_pad(2:end-1,2:end-1);
%bordermask=imdilate(bordermask,strel('disk',1,8)); %required because bwlabel can't discriminate a one-pixel separation
%bordermask=bwmorph(bordermask,'diag');
bordermask=bordermask(2:end-1,2:end-1);
DAs_pad=DAs_pad & ~bordermask;
DAs_pad=~bwmorph(~DAs_pad,'diag');
%DAs_pad=~bwmorph(~DAs_pad,'diag'); %clear remaining strands
end