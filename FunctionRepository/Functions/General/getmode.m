function modecalc=getmode(vals,numbins)
bmin=min(vals); bmax=max(vals); bstep=(bmax-bmin)/numbins; bins=bmin:bstep:bmax;
[kval,xval]=ksdensity(vals,bins);
maxidx=find(kval==max(kval),1); %first mode
modecalc=xval(maxidx);
end