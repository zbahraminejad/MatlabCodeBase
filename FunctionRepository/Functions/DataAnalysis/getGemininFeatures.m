function [risetime,badtraces]=getCdk2features(sampletraces,samplestats,minlength,relative)
slopewindow=min([minlength 40]); %default 40
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    signal=signal(firstframe(i):lastframe(i));
    %signal=smooth(signal(firstframe(i):lastframe(i)))';
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %sig_fwdslope=getslope_forward_avg(signal,round(slopewindow/2):slopewindow);
    sig_fwdslope=getslope_forward_avg(signal,6:10);
    sig_height=signal;
    heightgate=sig_height<0.3;
    filter=sig_fwdslope-sig_height+10;
    filter=filter.*heightgate;
    filter(1)=min(filter); %remove noisy first signal
    risetime(i)=find(filter==max(filter),1,'last');
    if filter(risetime(i))==0
        badtraces(i)=1;
    end
    %%% collect additional info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if relative==0
        risetime(i)=risetime(i)+firstframe(i)-1;
    end
    altstore1(i,firstframe(i):lastframe(i))=sig_fwdslope;
end
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%for i=1:samplesize
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    if badtraces(i)==1 || isnan(risetime(i))
        continue;
    end
    hold on;
    %plot(1:length(frames),altstore1(i,frames),'r');
    plot(risetime(i),sigstore(i,firstframe(i)+risetime(i)-1),'go','markerfacecolor','g','markersize',6);
end
%}