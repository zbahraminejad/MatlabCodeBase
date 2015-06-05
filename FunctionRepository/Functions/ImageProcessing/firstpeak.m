function ksthresh = firstpeak(values)
kmin=min(values); kmax=max(values); kstep=(kmax-kmin)/100;
[ksn,ksx]=ksdensity(values,kmin:kstep:kmax);
ksthresh=find(diff(ksn)<0,1,'first');
ksthresh=(ksthresh-1)/100;
