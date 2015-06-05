function [risetime,badtraces]=getPipdegRise_2(sampletraces,samplestats,minlength)
slopewindow=5;
gatewindow=min([minlength 5]);
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
risetimedb=risetime;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=sigstore;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    numframes=lastframe(i)-firstframe(i)+1;
    %signal=signal(firstframe(i):lastframe(i));
    signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% general requirements
    gate=zeros(1,numframes);
    for j=2:numframes-gatewindow
        pastheight=min(signal(1:j))<0.2 & max(signal(1:j))>0.3;
        futureheight=min(signal(j+1:end))>signal(j) & max(signal(j+gatewindow-1:end))>signal(j)+0.1;
        gate(j)=pastheight & futureheight;
    end
    sig_height=signal;
    sig_time=(1:length(signal))/length(signal);
    sig_revslope=getslope_reverse(signal,1:slopewindow);
    sig_fwdslope=getslope_forward(signal,1:slopewindow);
    %filter=sig_fwdslope-abs(sig_revslope)-2*sig_height+sig_time;
    filter=sig_fwdslope-abs(sig_revslope)-3*sig_height+sig_time;
    filter=filter.*gate;
    degend=find(filter(slopewindow+1:end)==max(filter(slopewindow+1:end)),1,'first')+slopewindow;
    altstore(i,samplestats(i,1):samplestats(i,2))=sig_fwdslope;
    %if sig_fwdslope(degend)<0.075 || sig_height(degend)>0.2
    %if sig_fwdslope(degend)<0.025 || sig_height(degend)>0.1
    deglast=min([degend+3 length(filter)]);
    %if max(sig_fwdslope(degend:deglast))<0.025 || sig_height(degend)>0.2
    if filter(degend)==0 || sig_height(degend)>0.2
        continue; %store as NaN to indicate that S-phase continues
    end
    risetimedb(i)=degend;
    risetime(i)=firstframe(i)+degend-1;
end
%keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    ylim([-0.1 2]);
    hold on;
    %plot(1:length(frames),altstore(i,frames),'r');
    if badtraces(i)==1
        framemid=round(median(frames));
        plot(framemid-frames(1)+1,sigstore(i,framemid),'rx','markersize',20);
        continue;
    elseif isnan(risetime(i))
        continue;
    end
    plot(risetimedb(i),sigstore(i,frames(risetimedb(i))),'go','markerfacecolor','g','markersize',6);
end
%}