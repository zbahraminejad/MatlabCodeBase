function [wellsss,blurframes]=removeblurs(wellsss,nucsizslot,fmax)

mednucsize=zeros(fmax,1);
for i=1:fmax
    mednucsize(i)=median(wellsss{i}(:,nucsizslot));
end
blurthresh=prctile(mednucsize,90)+3*iqr(mednucsize);
blurframes=find(mednucsize>blurthresh);
wellsss(blurframes)=[]; %if blurframes is empty, nothing happens
if isempty(blurframes)
    blurframes=0;
end

end