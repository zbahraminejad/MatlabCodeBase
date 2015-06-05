function [samplecells,samplestats]=gatetraces_cdk2(samplecells,DHBnuc,DHBratio,tracestats,window,threshold)
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%longtraces=tracestats(:,3)-tracestats(:,1)>=window;
longtraces=tracestats(:,3)>window; %default >=
samplecells=samplecells(longtraces);
samplestats=tracestats(longtraces,:);
DHBnuc=DHBnuc(longtraces,:);
DHBratio=DHBratio(longtraces,:);
%%% gate traces with low max DHBnuc intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
positivevals=min(DHBnuc,[],2)>threshold;
%positivevals=nanmedian(DHBnuc,2)>threshold;
samplecells=samplecells(positivevals);
samplestats=samplestats(positivevals,:);
DHBnuc=DHBnuc(positivevals,:);
DHBratio=DHBratio(positivevals,:);
%samplesignals(samplesignals<1)=1;
%%% gate DHB ratio outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    %tracenan=isnan(samplesignals2(i,:));
    %smoothtrace=smooth(samplesignals2(i,:));
    %smoothtrace(tracenan)=NaN;
    %samplesignals2(i,:)=smoothtrace;
    realframes=find(~isnan(DHBratio(i,:)));
    DHBratio(i,realframes)=smooth(DHBratio(i,realframes));
end
sig2min=min(DHBratio,[],2);
sig2max=max(DHBratio,[],2);
mu1=nanmean(DHBratio(:)); std1=nanstd(DHBratio(:));
iqr1=prctile(DHBratio(:),75)-prctile(DHBratio(:),25);
%outliers=sig2min<(mu1-3*std1) | sig2max>(mu1+3*std1);
outliers=sig2min<mu1-3*iqr1 | sig2max>mu1+3*iqr1;
%outliers=sig2min<0 | sig2max>5;
samplecells=samplecells(~outliers);
samplestats=samplestats(~outliers,:);
%DHBnuc=DHBnuc(~outliers,:);
DHBratio=DHBratio(~outliers,:);
%%% gate noisy DHBratio traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    %maxval=max(DHBnuc(i,:));
    %DHBnuc(i,:)=DHBnuc(i,:)/maxval;
    maxval=max(DHBratio(i,:));
    DHBratio(i,:)=DHBratio(i,:)/maxval;
end
swingthresh=1; %default 1
%noise=max(abs(diff(DHBnuc,2,2)),[],2)>swingthresh;
noise=max(abs(diff(DHBratio,2,2)),[],2)>swingthresh;
samplecells=samplecells(~noise);
samplestats=samplestats(~noise,:);