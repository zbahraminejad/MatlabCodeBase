function samplecells=gatetraces_cdk2(samplecells,samplesignals,tracestats,window)
%%% gate traces that are too short %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
longtraces=tracestats(:,3)-tracestats(:,1)>=window;
samplecells=samplecells(longtraces);
samplesignals=samplesignals(longtraces,:);
%%% gate traces without positive values and set floor to one %%%%%%%%%%%%%%
positivevals=max(samplesignals,[],2)>10;
samplecells=samplecells(positivevals);
samplesignals=samplesignals(positivevals,:);
samplesignals(samplesignals<1)=1;
%%% gate noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    maxval=max(samplesignals(i,:));
    samplesignals(i,:)=samplesignals(i,:)/maxval;
end
noise=max(abs(diff(samplesignals,2,2)),[],2)>1;
samplecells=samplecells(~noise);