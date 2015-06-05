function [tracesCdk2,badtracesCdk2]=gate_Cdk2_4(tracedata,nucchannel,cytochannel,tracestats,minlength,maxoption,maxthresh,noisethresh,quiescentanalysis)
numtraces=size(tracedata,1);
tracesCdk2gating=tracedata(:,:,nucchannel);
tracesCdk2=tracedata(:,:,cytochannel)./tracedata(:,:,nucchannel);
for i=1:numtraces
    realframes=find(~isnan(tracesCdk2gating(i,:)));
    tracesCdk2gating(i,realframes)=smooth(tracesCdk2gating(i,realframes));
    tracesCdk2(i,realframes)=smooth(tracesCdk2(i,realframes));
end
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtracesCdk2=gate_Cdk2_noisecalcdiff1_5(tracesCdk2gating,tracesCdk2,tracestats,minlength,maxoption,maxthresh,noisethresh,quiescentanalysis); %noisethresh=0.15
tracesCdk2=tracedata(:,:,cytochannel)./tracedata(:,:,nucchannel); %return raw signal (not smoothened)

