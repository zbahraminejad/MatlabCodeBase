function [tracesCdk2,badtracesCdk2]=gate_Cdk2(tracedata,nucchannel,cytochannel,tracestats,minlength,maxthresh,noisethresh)
numtraces=size(tracedata,1);
tracesCdk2gating=tracedata(:,:,nucchannel);
tracesCdk2=tracedata(:,:,cytochannel)./tracedata(:,:,nucchannel);
for i=1:numtraces
    realframes=find(~isnan(tracesCdk2gating(i,:)));
    tracesCdk2gating(i,realframes)=smooth(tracesCdk2gating(i,realframes));
    tracesCdk2(i,realframes)=smooth(tracesCdk2(i,realframes));
end
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtracesCdk2=gate_Cdk2_noisecalcdiff2(DHBnuc,tracesCdk2,daughterstats,minlengthdaughter,noisethresh); %noisethresh=1
badtracesCdk2=gate_Cdk2_noisecalcdiff1(tracesCdk2gating,tracesCdk2,tracestats,minlength,maxthresh,noisethresh); %noisethresh=0.15
