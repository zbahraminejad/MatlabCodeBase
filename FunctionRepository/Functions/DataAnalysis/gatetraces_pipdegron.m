function [samplecells,samplesignals,tracestats]=gatetraces_pipdegron(samplecells,samplesignals,tracestats)
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
longtraces=tracestats(:,3)-tracestats(:,1)>10;
samplecells=samplecells(longtraces);
samplesignals=samplesignals(longtraces,:);
tracestats=tracestats(longtraces,:);
%%% gate traces without positive values and set floor to zero %%%%%%%%%%%%%
positivevals=max(samplesignals,[],2)>50; %default 0. tried 50.
samplecells=samplecells(positivevals);
samplesignals=samplesignals(positivevals,:);
tracestats=tracestats(positivevals,:);
samplesignals(samplesignals<0)=0;
%%% normalize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    maxval=max(samplesignals(i,:));
    samplesignals(i,:)=samplesignals(i,:)/maxval;
end
%%% gate noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise=max(abs(diff(samplesignals,2,2)),[],2)>1;
samplecells=samplecells(~noise);
samplesignals=samplesignals(~noise,:);
tracestats=tracestats(~noise,:);