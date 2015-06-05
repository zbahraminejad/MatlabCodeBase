function [degstarts,degends,badtraces]=getdegstartandend_openended(sampletraces,relative)
samplesize=size(sampletraces,1);
degstarts=ones(samplesize,1)*NaN;
degends=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
%sigstore=ones(samplesize,tracelength)*NaN;
%altstore=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe=find(~isnan(signal),1,'first');
    lastframe=find(~isnan(signal),1,'last');
    %signal=smooth(signal(firstframe:lastframe))';  %default was to smooth this!
    signal=signal(firstframe:lastframe);
    %sigstore(i,:)=signal;
    %%% calc degstart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig_revslope=getslope_reverse(signal,1:5); %default 10
    sig_height=signal;
    filter=sig_revslope-sig_height;
    lowerleft=find(filter==max(filter),1,'first');
    if lowerleft<=7  %noisy signal
        badtraces(i)=1;
        continue;
    elseif filter(lowerleft)<0.25  %hasn't entered S-phase yet
        continue;
    end
    sig_time=(1:length(signal))/length(signal);
    sig_fwdslope=getslope_forward(signal,1:3);
    filter=sig_time-sig_revslope-sig_fwdslope*2;
    %degstarts(i)=find(filter(1:lowerleft-1)==max(filter(1:lowerleft-1)),1,'last')+1;
    degstarts(i)=lowerleft;
    if degstarts(i)<5
        badtraces(i)=1;
    end
    if relative==0
        degstarts(i)=degstarts(i)+firstframe-1;
    end;
    %%% calc degend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %signal=smooth(signal)';
    sig_revslope=getslope_reverse(signal,1:5);
    sig_fwdslope=getslope_forward(signal,1:5);
    filter=sig_fwdslope-abs(sig_revslope)-2*sig_height+sig_time;
    lowerright=find(filter(lowerleft+1:end)==max(filter(lowerleft+1:end)))+lowerleft;
    if lowerleft==length(signal) || lowerright==length(signal)
        continue;
    elseif sig_fwdslope(lowerright)<0.01
        badtraces(i)=1;
        continue;
    end
    if relative==0
        degends(i)=lowerright+firstframe-1;
    else
        degends(i)=lowerright;
    end
    %altstore(i,:)=sig_fwdslope;
end
badtraces=badtraces>0;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
x=1:length(signal);
for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    plot(x,sigstore(i,:));
    xlim([1 length(signal)]);
    if badtraces(i)==1
        continue;
    end
    hold on;
    %plot(x,altstore(i,:),'r');
    plot(degstarts(i),sigstore(i,degstarts(i)),'go','markerfacecolor','g','markersize',6);
    plot(degends(i),sigstore(i,degends(i)),'ro','markerfacecolor','r','markersize',6);
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