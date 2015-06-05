function [samplecells,samplestats]=gatetraces_cdk2(samplecells,samplesignals,tracestats,window,threshold)
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%longtraces=tracestats(:,3)-tracestats(:,1)>=window;
longtraces=tracestats(:,3)>window; %default >=
samplecells=samplecells(longtraces);
samplestats=tracestats(longtraces,:);
samplesignals=samplesignals(longtraces,:);
%%% gate traces without positive values and set floor to one %%%%%%%%%%%%%%
positivevals=min(samplesignals,[],2)>threshold;
%positivevals=nanmedian(samplesignals,2)>threshold;
samplecells=samplecells(positivevals);
samplestats=samplestats(positivevals,:);
samplesignals=samplesignals(positivevals,:);
%samplesignals(samplesignals<1)=1;
%%% gate ratio outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smoothsamplesignals=ones(size(samplesignals))*NaN;
for i=1:numel(samplecells)
    %tracenan=isnan(samplesignals(i,:));
    %smoothtrace=smooth(samplesignals(i,:));
    %smoothtrace(tracenan)=NaN;
    %smoothsamplesignals(i,:)=smoothtrace;
    realframes=find(~isnan(samplesignals(i,:)));
    smoothsamplesignals(i,realframes)=smooth(samplesignals(i,realframes));
end
sigmin=min(smoothsamplesignals,[],2);
sigmax=max(smoothsamplesignals,[],2);
mu1=nanmean(smoothsamplesignals(:)); std1=nanstd(smoothsamplesignals(:));
outliers=sigmin<(mu1-3*std1) | sigmax>(mu1+3*std1);
samplecells=samplecells(~outliers);
samplestats=samplestats(~outliers,:);
samplesignals=samplesignals(~outliers,:);
%%% gate noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    maxval=max(samplesignals(i,:));
    samplesignals(i,:)=samplesignals(i,:)/maxval;
end
swingthresh=1; %default 1
noise=max(abs(diff(samplesignals,2,2)),[],2)>swingthresh;
samplecells=samplecells(~noise);
samplestats=samplestats(~noise,:);