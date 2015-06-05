function [degstarts,badtraces]=getdegstarts(sampletraces)
[samplesize,tracelength]=size(sampletraces);
degstarts=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    sigstore(i,:)=signal;
    sig_revslope=getslope_reverse(signal,1:10);
    sig_height=-signal;
    filter=sig_revslope+sig_height;
    lowerleft=find(filter==max(filter),1,'first');
    if filter(lowerleft)<0.25
        badtraces(i)=1;
        degstarts(i)=1;
        continue;
    end
    sig_time=(1:length(signal))/length(signal);
    sig_fwdslope=getslope_forward(signal,1:3);
    filter=sig_time-sig_revslope-sig_fwdslope*2;
    degstarts(i)=find(filter(1:lowerleft-1)==max(filter(1:lowerleft-1)),1,'last');
    altstore(i,:)=filter;
end
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
x=1:length(signal);
for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    plot(x,sigstore(i,:));
    hold on;
    %plot(x,altstore(i,:),'r');
    plot(degstarts(i),sigstore(i,degstarts(i)),'ro','markerfacecolor','r','markersize',6);
    xlim([1 length(signal)]);
end
%}