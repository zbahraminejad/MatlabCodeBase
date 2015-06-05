function [tracesCdk2,badtracesCdk2]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis,removequiescent)
numtraces=size(tracedata,1);
DHBnuc=tracedata(:,:,nucchannel);
DHBratio=tracedata(:,:,cytochannel)./tracedata(:,:,nucchannel);
for i=1:numtraces
    realframes=find(~isnan(DHBnuc(i,:)));
    DHBnuc(i,realframes)=smooth(DHBnuc(i,realframes));
    DHBratio(i,realframes)=smooth(DHBratio(i,realframes));
end
%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracestats,1);
if nucthreshoption==0
%     lowsignal=max(DHBnuc,[],2)<threshold;
    lowsignal=prctile(DHBnuc,75,2)<nucthresh;
elseif nucthreshoption==1
    firstval=ones(numtraces,1)*NaN;
    for i=1:numtraces
        firstval(i)=DHBnuc(i,tracestats(i,1));
    end
    lowsignal=firstval<nucthresh;
end
%%% remove ratio outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig2min=min(DHBratio,[],2);
sig2max=max(DHBratio,[],2);
%mu1=nanmean(DHBratio(:));
mumin=nanmean(sig2min);
mumax=nanmean(sig2max);
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
%%% combine all gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if quiescentanalysis
    badtracesCdk2=lowsignal | outliers | noisy;
    notquiescent=mean(DHBratio(:,1:3),2)>0.4 | tracestats(:,1)>1 | min(DHBratio,[],2)>0.5 | max(DHBratio,[],2)>0.5; %only accept cells that are present from beginning and <1.0 ratio
    badtracesCdk2=badtracesCdk2 | notquiescent;
elseif removequiescent
    cyclingcells = max(DHBratio,[],2)>1;
    badtracesCdk2 = lowsignal | outliers | noisy | ~cyclingcells;
    
elseif motherthresh>0
    % remove traces where max DHB ratio in mother is too low
    mothermaxval=ones(numtraces,1)*NaN;
    for i=1:numtraces
        mothermaxval(i)=max(DHBratio(i,1:tracestats(i,1)-1));
    end
    lowmother=mothermaxval<motherthresh;
    badtracesCdk2=lowsignal | outliers | noisy | lowmother;
else
    badtracesCdk2=lowsignal | outliers | noisy;
end

tracesCdk2=tracedata(:,:,cytochannel)./tracedata(:,:,nucchannel); %return raw signal (not smoothened)


%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(max(DHBnuc,[],2),0:40:4040); xlim([0 4000]);
%hist(prctile(DHBnuc,75,2),-20:20:1020); xlim([0 1000]);
%hist(firstval,2000);
%hist(noise,-2.01:0.01:2.01); xlim([-2 2]);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%