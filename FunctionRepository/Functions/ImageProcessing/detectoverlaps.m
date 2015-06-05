function excludemap = detectoverlaps(home,out,orientation,cytoringoutermass,legitedges,imap,height)
pixelidx=regionprops(home,'PixelIdxList');
crossnuc=out.*cytoringoutermass;
crossgrad=out.*legitedges;
uniquecn=unique(crossnuc);uniquecn=uniquecn(uniquecn>0);
uniquecg=unique(crossgrad);uniquecg=uniquecg(uniquecg>0);
crossboth=ismember(uniquecn,uniquecg);
overlap=uniquecn(~crossboth);
potentialoverlap=uniquecn(crossboth);
verifiedoverlap=zeros(size(potentialoverlap));
crossnuc=regionprops(crossnuc,'PixelIdxList'); %returns cell ID of any spoke hitting neighbor nucleus
for k=1:length(crossnuc)
    crossnuc(k).PixelIdxList=abs(max(orientation*crossnuc(k).PixelIdxList));
end
crossgrad=regionprops(crossgrad,'PixelIdxList');  %returns cell ID of any spoke hitting cyto gradient
for k=1:length(crossgrad)
    crossgrad(k).PixelIdxList=abs(max(orientation*crossgrad(k).PixelIdxList));
end
for k=1:length(potentialoverlap)
    kval=potentialoverlap(k);
    if orientation*crossgrad(kval).PixelIdxList<orientation*crossnuc(kval).PixelIdxList  %gradient cross is north of nucleus cross
        verifiedoverlap(k)=1;
    end
end
overlap=[overlap;potentialoverlap(logical(verifiedoverlap))];
trueoverlap=zeros(size(overlap));
for k=1:length(overlap)
    kval=overlap(k);
    neix=crossnuc(kval).PixelIdxList;
    neic=[mod(neix,height) ceil(neix/height)];
    nucx=pixelidx(kval).PixelIdxList;
    nucc=[mod(nucx,height) ceil(nucx/height)];
    bothc=[neic;nucc];
    bothc(bothc==0)=height;
    midc=round(mean(bothc,1));
    if imap(neix)>imap(nucx) && imap(midc(1),midc(2))>imap(nucx)    %both the neighbor and midpoint are higher than ring
        trueoverlap(k)=1;
    end
end
overlap=overlap(logical(trueoverlap));
excludemap=ismember(home,overlap);
end