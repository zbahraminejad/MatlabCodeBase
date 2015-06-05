function badtraces=gate_Cdk2(DHBnuc,DHBratio,tracestats,minlength,threshold,noisethresh,quiescentanalysis)
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(max(DHBnuc,[],2),0:40:4040); xlim([0 4000]);
%hist(prctile(DHBnuc,75,2),0:8:1008); xlim([0 1000]);
%hist(noise,-2.01:0.01:2.01); xlim([-2 2]);
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<minlength; %as of 2014-03-17 was <=
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
lowsignal=max(DHBnuc,[],2)<threshold; %as of 2014-03-17 was <=
%%% remove ratio outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig2min=min(DHBratio,[],2);
sig2max=max(DHBratio,[],2);
%mu1=nanmean(DHBratio(:));
mumin=mean(sig2min);
mumax=mean(sig2max);
%iqr1=prctile(DHBratio(:),75)-prctile(DHBratio(:),25);
%outliers=sig2min<mu1-3*iqr1 | sig2max>mu1+3*iqr1;
iqrmin=prctile(sig2min,75)-prctile(sig2min,25);
iqrmax=prctile(sig2max,75)-prctile(sig2max,25);
outliers=sig2min<mumin-3*iqrmin | sig2max>mumax+3*iqrmax;
%%% remove noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%noisethreshthresh=0.15; %default 0.15
%maxidxmult=floor((numframes-1)/5);
%sampleindices=1:5:(maxidxmult*5+1);
%shortDHBratio=DHBratio(:,sampleindices);
noise=max(diff(DHBratio,1,2),[],2);
noisy=noise>noisethresh;

badtraces=shorttraces | lowsignal | outliers | noisy;
if quiescentanalysis
    notquiescent=mean(DHBratio(:,1:3),2)>1.0 | tracestats(:,1)>1; %only accept cells that are present from beginning and <1.0 ratio
    badtraces=badtraces | notquiescent;
end

