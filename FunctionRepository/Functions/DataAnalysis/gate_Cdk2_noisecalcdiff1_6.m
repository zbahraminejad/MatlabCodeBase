function badtraces=gate_Cdk2(DHBnuc,DHBratio,tracestats,minlength,maxoption,threshold,maxmother,noisethresh,quiescentanalysis)
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(max(DHBnuc,[],2),0:40:4040); xlim([0 4000]);
%hist(prctile(DHBnuc,75,2),-20:20:1020); xlim([0 1000]);
%hist(firstval,2000);
%hist(noise,-2.01:0.01:2.01); xlim([-2 2]);
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shorttraces=tracestats(:,3)<minlength; %as of 2014-03-17 was <=
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracestats,1);
if maxoption==0
    %lowsignal=max(DHBnuc,[],2)<threshold;
    lowsignal=prctile(DHBnuc,75,2)<threshold;
elseif maxoption==1
    firstval=ones(numtraces,1)*NaN;
    for i=1:numtraces
        firstval(i)=DHBnuc(i,tracestats(i,1));
    end
    lowsignal=firstval<threshold;
end
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
noise=ones(numtraces,1)*NaN;
for i=1:numtraces
    maxdiff=max(diff(DHBratio(i,tracestats(i,1):tracestats(i,2)),1));
    if ~isempty(maxdiff)
        noise(i)=maxdiff;
    end
end
noisy=noise>noisethresh;

if quiescentanalysis
    badtraces=shorttraces | lowsignal | outliers | noisy;
    %notquiescent=mean(DHBratio(:,1:3),2)>0.75 | tracestats(:,1)>1; %only accept cells that are present from beginning and <1.0 ratio
    notquiescent=mean(DHBratio(:,1:3),2)>0.75 | tracestats(:,1)>1 | min(DHBratio,[],2)>0.5; %only accept cells that are present from beginning and <1.0 ratio
    badtraces=badtraces | notquiescent;
else
    %%% remove traces where max DHB ratio in mother is too low %%%%%%%%%%%%
    mothermaxval=ones(numtraces,1)*NaN;
    for i=1:numtraces
        mothermaxval(i)=max(DHBratio(i,1:tracestats(i,1)-1));
    end
    lowmother=mothermaxval<maxmother;
    badtraces=shorttraces | lowsignal | outliers | noisy | lowmother;
end

