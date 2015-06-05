function badtraces=gate_Cdk2(DHBnuc,DHBratio,tracestats,minlength,threshold,noisethresh)
numtraces=size(tracestats,1);
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<=minlength;
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
lowsignal=min(DHBnuc,[],2)<=threshold;
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
for i=1:numtraces
    maxval=max(DHBratio(i,:));
    DHBratio(i,:)=DHBratio(i,:)/maxval;
end
noisy=max(abs(diff(DHBratio,2,2)),[],2)>noisethresh;

badtraces=shorttraces | lowsignal | outliers | noisy;