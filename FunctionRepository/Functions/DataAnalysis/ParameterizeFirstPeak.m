function [peakstart,lowerthresh,peakmid,upperthresh,peakend,lastval]=ParameterizeFirstPeak(vector,subprctile)
vector=vector(:);
lowerbound=min(vector);
%upperbound=prctile(vector,99);
upperbound=1500000;
%bin=(upperbound-lowerbound)/50;
%[kn,kx]=ksdensity(vector,lowerbound+bin/2:bin:upperbound-bin/2);
bstep=(upperbound-lowerbound)/100;
bin=lowerbound:bstep:upperbound;
n_elements = histc(vector,bin);
n_elements= 100*n_elements/sum(n_elements);
histdiff=diff(n_elements);
peakstart=find(histdiff>0.5,1,'first');
histrise=histdiff>=0;
histrise(1:peakstart)=1;
peakmid=find(histdiff<-1,1,'first');
histrise(1:peakmid)=0;
%peakend=find(histrise,1,'first');
zerovals=n_elements(1:end-1)==0 & histrise;
peakend=find(zerovals,1,'first');

countedvals=n_elements>0;
countedvals=imerode(countedvals,[1;1;1]);
lastval=find(countedvals,1,'last');

peakcdf=cumsum(n_elements(peakstart:peakend));
peakend=peakstart-1+find(peakcdf>=90,1,'first');
lowerthresh=peakstart-1+find(peakcdf>=subprctile,1,'first');
upperthresh=peakstart-1+find(peakcdf>=100-subprctile,1,'first');

%%% convert to vector vals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
peakstart=bin(peakstart);
peakmid=bin(peakmid);
peakend=bin(peakend);
lowerthresh=bin(lowerthresh);
upperthresh=bin(upperthresh);
lastval=bin(lastval);