function [samplecells,samplesignals]=gatetraces_geminin(samplecells,samplesignals,samplestats,sigmax,window)
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
longtraces=samplestats(:,3)-samplestats(:,1)>=window;
samplecells=samplecells(longtraces);
samplesignals=samplesignals(longtraces,:);
sigmax=sigmax(longtraces);
%%% gate traces without positive values and set floor to one %%%%%%%%%%%%%%
daughtermax=max(samplesignals,[],2);
sigmax=max([sigmax daughtermax],[],2);
positivevals=sigmax>0;
samplecells=samplecells(positivevals);
samplesignals=samplesignals(positivevals,:);
samplesignals(samplesignals<1)=1;
sigmax=sigmax(positivevals);
%%% normalize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    samplesignals(i,:)=samplesignals(i,:)/sigmax(i);
end
%%% gate noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise=max(abs(diff(samplesignals,2,2)),[],2)>1;
samplecells=samplecells(~noise);
samplesignals=samplesignals(~noise,:);