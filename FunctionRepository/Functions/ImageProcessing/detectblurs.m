function blurframes=detectblurs(wellsss,nucsizslot)
fmax=size(wellsss,3);
mednucsize=zeros(fmax,1);
for i=1:fmax
    mednucsize(i)=median(wellsss{i}(:,nucsizslot));
end
blurthresh=prctile(mednucsize,90)+3*iqr(mednucsize);
blurframes=find(mednucsize>blurthresh);
if isempty(blurframes)
    blurframes=0;
end

end