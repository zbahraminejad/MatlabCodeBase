function [falltime,badtraces]=getPipdegFall(sampletraces,samplestats,minlength)
gatewindow=min([minlength 5]);
%This function will return degstart assuming trace begins with mitosis
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
    numframes=lastframe(i)-firstframe(i)+1;
    signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%%%%% general requirements
    gate=zeros(1,numframes); gate(1)=1;
    for j=2:numframes-gatewindow
        pastheight=min(signal(1:j))>0.1;
        futureheight=max(signal(j+1:j+gatewindow))<signal(j);
        gate(j)=pastheight && futureheight;
    end
    diffsig=diff(signal,1);
    dssmooth=smooth(diffsig)';
    gate=gate(1:end-1);
    dssmooth=dssmooth.*gate;
    if min(dssmooth)>-0.05
        if min(signal)<0.1
            badtraces(i)=1;
        end
        continue;
    end
    maxdropidx=find(dssmooth==min(dssmooth),1,'first');
    dssmoothneg=dssmooth<0;
    dssnlabel=bwlabel(dssmoothneg);
    maxdroplabel=dssnlabel(maxdropidx);
    degstart=find(dssnlabel==maxdroplabel,1,'first');
    falltimedb(i)=degstart;
    falltime(i)=firstframe(i)+degstart-1;
end
%keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames)); ylim([-0.5 3]);
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