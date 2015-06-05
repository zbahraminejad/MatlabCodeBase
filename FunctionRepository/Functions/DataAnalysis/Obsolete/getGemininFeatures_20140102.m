function [risetime,badtraces]=getCdk2features(sampletraces,samplestats,minlength,relative)
slopewindow=min([minlength 40]); %default 40
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
altvar1=ones(samplesize,1)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    signal=smooth(signal(firstframe(i):lastframe(i)))';
    numframes=lastframe(i)-firstframe(i)+1;
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    badtraces(i)=min(signal(1:slopewindow))>1; %will remove permanently high signals
    %%%%%% filter for inflection points
    sig_fwdslope_short=getslope_forward_avg(signal,6:10);
    sig_fwdslope_long=getslope_forward_avg(signal,round(slopewindow/2+1):slopewindow);
    sig_time=(1:length(signal))/length(signal);
    filter=sig_fwdslope_short+sig_fwdslope_long-2*signal-sig_time*0.01;
    %%% find cdk2start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filtermax=max(filter);
    risetime(i)=find(filter==filtermax,1,'first');
    if max(signal)<=0.2
        risetime(i)=NaN;
        continue;
    end
    %%% collect additional info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    altvar1(i)=sig_fwdslope_long(risetime(i));
    if relative==0
        risetime(i)=risetime(i)+firstframe(i)-1;
    end
    altstore1(i,firstframe(i):lastframe(i))=filter;
end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%for i=1:samplesize
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    %framelength=lastframe(i)-firstframe(i)+1;
    %frames=1:framelength;
    %plot(frames,sigstore(i,frames));
    %axis([frames(1) frames(end) 0 2.5]);
    plot(1:length(frames),sigstore(i,frames));
    %axis([1 length(frames) 0 1]);
    if badtraces(i)==1 || isnan(risetime(i))
        continue;
    end
    hold on;
    plot(1:length(frames),altstore1(i,frames),'r');
    %plot(1:length(frames),altstore1(i,frames),'r');
    plot(risetime(i),sigstore(i,risetime(i)),'go','markerfacecolor','g','markersize',6);
    %plot(risetime(i):risetime(i)+9,minval(i)+riseslope(i)*(0:9),'r');
end
%}