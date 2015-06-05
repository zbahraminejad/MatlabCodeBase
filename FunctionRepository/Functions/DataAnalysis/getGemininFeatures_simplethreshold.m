function [risetime,badtraces]=getCdk2features(sampletraces,samplestats,minlength,relative)
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    signal=signal(firstframe(i):lastframe(i));
    %signal=smooth(signal(firstframe(i):lastframe(i)))';
    sigstore(i,firstframe(i):lastframe(i))=signal;
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prebuffer=5;
    badtraces(i)=signal(prebuffer+1)>0.02;
    tempsearch=find(signal(prebuffer+1:end)>0.02,1,'first');
    if isempty(tempsearch)
        risetime(i)=NaN;
    else
        risetime(i)=prebuffer+tempsearch;
    end
end
keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:96
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(i):lastframe(i);
    plot(1:length(frames),sigstore(i,frames));
    if badtraces(i)==1 || isnan(risetime(i))
        continue;
    end
    hold on;
    plot(risetime(i),sigstore(i,frames(risetime(i))),'go','markerfacecolor','g','markersize',6);
end
%}