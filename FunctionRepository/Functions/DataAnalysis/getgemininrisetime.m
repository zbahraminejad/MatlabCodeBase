function [risetime,badtraces]=getgemininrisetime(sampletraces,relative)
[samplesize,tracelength]=size(sampletraces);
risetime=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
sigstore=ones(samplesize,tracelength)*NaN;
altstore=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=find(~isnan(signal),1,'first');
    lastframe(i)=find(~isnan(signal),1,'last');
    signal=smooth(signal(firstframe(i):lastframe(i)))';  %default was to smooth this!
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% calc geminin risetime %%%%%%%%%%%%%%%%%%%%%%%
    sig_fwdslope=getslope_forward_avg(signal,6:10);
    sig_height=signal;
    heightgate=sig_height<0.1;
    %filter=2*sig_fwdslope-sig_height+sig_time+2;
    filter=50*sig_fwdslope-5*sig_height+0.5;
    filter=filter.*heightgate;
    filter(1)=min(filter); %remove noisy first signal
    %filter((end-9):end)=min(filter);
    risetime(i)=find(filter==max(filter),1,'last');
    if min(signal(risetime(i):end))<signal(risetime(i)) || signal(end)<0.1
        badtraces(i)=1;
    end
    if relative==0
        risetime(i)=risetime(i)+firstframe(i)-1;
    end
    altstore(i,firstframe(i):lastframe(i))=filter;
end
badtraces=badtraces>0;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%for i=1:samplesize
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(frames,sigstore(i,frames));
    axis([frames(1) frames(end) 0 2]);
    if badtraces(i)==1
        continue;
    end
    hold on;
    plot(frames,altstore(i,frames),'r');

    %plot(degstarts(i),sigstore(i,degstarts(i)),'go','markerfacecolor','g','markersize',6);
    %plot(degends(i),sigstore(i,degends(i)),'ro','markerfacecolor','r','markersize',6);
end


plot(x,sigstore(22,:),'linewidth',4); xlim([1 length(signal)]);;
hold on;
plot(degstarts(22),sigstore(22,degstarts(22)),'go','markerfacecolor','g','markersize',12);
plot(degends(22),sigstore(22,degends(22)),'ro','markerfacecolor','r','markersize',12);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');

good=find(badtraces==0 & degends<length(signal));
bad=find(badtraces==0 & degends==length(signal));
goodvals=ones(numel(good),1)*NaN;
for i=1:numel(good)
    g=good(i);
    goodvals(i)=altstore(g,degends(g));
end
figure,hist(goodvals);
badvals=ones(numel(bad),1)*NaN;
for i=1:numel(bad)
    b=bad(i);
    badvals(i)=altstore(b,degends(b));
end
figure,hist(badvals);
%}