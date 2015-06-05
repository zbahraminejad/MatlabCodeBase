function [risetime,badtraces]=getCdk2features_steve(sampletraces,samplestats,minlength)
slopewindow=min([minlength 10]); %default 40
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
%minval=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
risetimedb=risetime;
%riseslope=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore1=ones(samplesize,tracelength)*NaN;
altvar1=ones(samplesize,1)*NaN;
for i=1:samplesize
    signal_total=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    %signal=signal(firstframe(i):lastframe(i));
    signal_smooth(i,:)=smooth(signal_total);
    signal=smooth(signal_total(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    numframes=lastframe(i)-firstframe(i)+1;
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     badtraces(i)=min(signal(1:slopewindow))>1; %will remove permanently high signals
    
    badtraces(i)=min(signal_total(firstframe(i):firstframe(i)+slopewindow))>1; %will remove permanently high signals
    
    
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
    filter=sig_fwdslope_short+sig_fwdslope_long-4*signal+sig_time+10;
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
    minval(i)=min(signal(1:risetime(i)));
    %%% calc riseslope %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %riseslope(i)=sig_fwdslope_short(risetime(i));
    %%% collect additional info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    altvar1(i)=sig_fwdslope_long(risetime(i));
    risetimedb(i)=risetime(i);
    risetime(i)=risetime(i)+firstframe(i); %return absolute POI rather than relative to mitosis
    altstore1(i,firstframe(i):lastframe(i))=filter;
    
%     if risetime(i)-firstframe(i)<length(signal)-20;
%         if signal_total(risetime(i)+20)-signal_total(risetime(i))<0.14;
%             risetime(i)=NaN;
%         end
%     end

end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
altstore1(altstore1==0)=NaN;
altstore_norm=normalizeMyTracesGeminin_alt2(altstore1,0.01);
for i=1:96
%for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    axis([1 length(frames) 0.3 2.5]);

    if badtraces(i)==1 || isnan(risetimedb(i))
        continue;
    end
    hold on;
    plot(risetimedb(i),sigstore(i,frames(risetimedb(i))),'go','markerfacecolor','g','markersize',6);
    %plot(1:length(frames),altstore_norm(i,frames)+0.5,'r');
    %plot(risetimedb(i),altstore_norm(i,frames(risetimedb(i)))+0.5,'go','markerfacecolor','g','markersize',6);
end
%}

%{
altstore1(altstore1==0)=NaN;
altstore_norm=normalizeMyTracesGeminin_alt2(altstore1,0.01);
for i=1:samplesize
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    numframes=size(altstore1,2);    
    frames=1:numframes;
    plot(1:length(frames),signal_smooth(i,frames));
    axis([1 length(frames) 0.3 2.5]);
    hold on;
    plot(1:length(frames),altstore_norm(i,frames)+0.5,'r');
    if badtraces(i)==1 || isnan(risetime(i))
        continue;
    end
    hold on;
    plot(risetime(i),signal_smooth(i,frames(risetime(i))),'go','markerfacecolor','g','markersize',6);
    plot(risetime(i),altstore_norm(i,frames(risetime(i)))+0.5,'go','markerfacecolor','g','markersize',6);
end
%}