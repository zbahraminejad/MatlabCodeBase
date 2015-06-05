function [minval,badtraces]=getCdk2exitval(sampletraces)
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
minval=ones(samplesize,1)*NaN;
mintime=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=find(~isnan(signal),1,'first');
    lastframe(i)=find(~isnan(signal),1,'last');
    signal=smooth(signal(firstframe(i):lastframe(i)))';
    heightgate=signal<1.1 & signal>0.2;
    if length(signal)>50
        heightgate(50:end)=0;
    end
    if max(heightgate)==0
        badtraces(i)=1;
        continue;
    end
    signaldepth=1.1-signal;
    signaldepth=signaldepth.*heightgate;
    signaldepth(1)=min(signaldepth); %avoid noisy first frame
    mintime(i)=find(signaldepth==max(signaldepth),1,'first');
    minval(i)=signal(mintime(i));
    sigstore(i,firstframe(i):lastframe(i))=signal;
    altstore1(i,firstframe(i):lastframe(i))=heightgate;
end
badtraces=badtraces>0;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%for i=1:samplesize
for i=97:192
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(frames,sigstore(i,frames));
    axis([frames(1) frames(end) 0 2]);
    if badtraces(i)==1
        continue;
    end
    hold on;
    %plot(frames,altstore1(i,frames),'r');
    %plot(firstframe(i)+mintime(i)-1,minval(i),'go','markerfacecolor','g','markersize',6);
    %plot(risetime(i):risetime(i)+9,minval(i)+riseslope(i)*(0:9),'r');
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