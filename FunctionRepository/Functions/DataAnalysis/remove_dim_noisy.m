function samplecells=remove_dim_noisy(samplesignals,samplecells,threshold)
%%% gate traces without positive values and set floor to one %%%%%%%%%%%%%%
%positivevals=max(samplesignals,[],2)>threshold;
positivevals=min(samplesignals,[],2)>threshold;
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