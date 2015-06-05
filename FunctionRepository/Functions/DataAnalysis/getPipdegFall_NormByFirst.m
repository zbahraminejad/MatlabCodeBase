function [falltime,badtraces]=getPipdegFall(sampletraces,samplestats,minlength)
%This function will return degstart assuming trace begins with mitosis
slopewindow=min([minlength 10]);
slopewindowshort=min([minlength 5]);
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
falltime=ones(samplesize,1)*NaN;
falltimedb=falltime;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=sigstore;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build falltime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig_height=signal;
    sig_time=(1:length(signal))/length(signal);
    sig_revslope=getslope_reverse(signal,1:slopewindow);
    sig_fwdslope=getslope_forward(signal,1:slopewindow);
    sig_revslopeshort=getslope_reverse(signal,1:slopewindowshort);
    filter=sig_revslopeshort-abs(sig_fwdslope)-2*sig_height-sig_time;
    degstart=find(filter(slopewindowshort+1:end)==max(filter(slopewindowshort+1:end)),1,'first')+slopewindowshort;
    altstore(i,firstframe(i):lastframe(i))=sig_revslope;
    if sum(sig_height(1:5)<0.1)>0 || abs(sig_height(degstart))>0.2
        badtraces(i)=1;
        continue;
    end
    if sig_revslope(degstart)<0.1
        continue;
    end
    falltimedb(i)=degstart;
    falltime(i)=firstframe(i)+degstart-1;
end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames)); ylim([0 1.5]);
    hold on;
    if badtraces(i)==1
        framemid=round(median(frames));
        plot(framemid-frames(1)+1,sigstore(i,framemid),'rx','markersize',20);
        continue;
    elseif isnan(falltime(i))
        continue;
    end
    plot(falltimedb(i),sigstore(i,frames(falltimedb(i))),'go','markerfacecolor','g','markersize',6);
    %plot(1:length(frames),altstore(i,frames),'r');
end
%}