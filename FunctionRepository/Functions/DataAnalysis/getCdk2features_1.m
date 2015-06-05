function [risetime,badtraces]=getCdk2features_1(sampletraces,samplestats,minlength)
slopewindow=min([minlength 40]); %default 40
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
%minval=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
%riseslope=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
altvar1=ones(samplesize,1)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    %signal=signal(firstframe(i):lastframe(i));
    signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    numframes=lastframe(i)-firstframe(i)+1;
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    badtraces(i)=min(signal(1:slopewindow))>1; %will remove permanently high signals
    %%%%%% general requirements
    gate=zeros(1,numframes);
    for j=2:numframes-slopewindow
        pastheight=max(signal(1:j))<1.5;
        futureheight=max(signal(j:end))>1;
        gate(j)=pastheight && futureheight;
    end
    %%%%%% filter for inflection points
    sig_fwdslope_short=getslope_forward_avg(signal,6:10);
    sig_fwdslope_long=getslope_forward_avg(signal,round(slopewindow/2+1):slopewindow);
    sig_time=(1:length(signal))/length(signal);
    filter=sig_fwdslope_short+sig_fwdslope_long-2*signal+sig_time+10;
    %%%%%% combine gate and filter
    filter=filter.*gate;
    %%% find cdk2start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filtermax=max(filter);
    risetime(i)=find(filter==filtermax,1,'first');
    if filtermax<=0
        risetime(i)=NaN;
        continue;
    end
    %%% calc minval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %minval(i)=min(signal(1:risetime(i)));
    %%% calc riseslope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %riseslope(i)=sig_fwdslope_short(risetime(i));
    %%% collect additional info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    altvar1(i)=sig_fwdslope_long(risetime(i));
    risetime(i)=risetime(i)+firstframe(i)-1; %return absolute POI rather than relative to mitosis
    altstore1(i,firstframe(i):lastframe(i))=filter;
end

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%for i=1:samplesize
for i=101:196
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    %plot(frames,sigstore(i,frames));
    %axis([frames(1) frames(end) 0 2.5]);
    plot(1:length(frames),sigstore(i,frames));
    axis([1 length(frames) 0 2.5]);
    if badtraces(i)==1
        continue;
    end
    hold on;
    %plot(frames,altstore1(i,frames),'r');
    %plot(1:length(frames),altstore1(i,frames),'r');
    plot(risetime(i),minval(i),'go','markerfacecolor','g','markersize',6);
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