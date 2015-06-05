function degends=getdegends(sampletraces)
[samplesize,tracelength]=size(sampletraces);
degends=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=ones(samplesize,tracelength)*NaN;
altstore2=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    sigstore(i,:)=signal;
    signal=smooth(signal)';
    sig_revslope=getslope_reverse(signal,1:5);
    sig_fwdslope=getslope_forward(signal,1:5);
    sig_height=signal;
    sig_time=(1:length(signal))/length(signal);
    filter=sig_fwdslope-abs(sig_revslope)-2*sig_height+sig_time;
    lowerright=find(filter==max(filter));
    degends(i)=lowerright;
    altstore(i,:)=filter;
    altstore2(i,:)=2*sig_fwdslope-abs(sig_revslope);
end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
x=1:length(signal);
for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    plot(x,sigstore(i,:));
    hold on;
    %plot(x,altstore(i,:),'r');
    %plot(x,altstore2(i,:),'g');
    plot(degends(i),sigstore(i,degends(i)),'ro','markerfacecolor','r','markersize',6);
    xlim([1 length(signal)]); ylim([-0.2 1]);
end
%}