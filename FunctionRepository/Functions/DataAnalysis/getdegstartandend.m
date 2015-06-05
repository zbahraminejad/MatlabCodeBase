function [degstarts,degends]=getdegstartandend_mother(sampletraces,stats,minlength,relative)
slopewindow=min([minlength 10]);
slopewindowshort=min([minlength 5]);
[samplesize,tracelength]=size(sampletraces);
degstarts=ones(samplesize,1)*NaN;
degends=ones(samplesize,1)*NaN;
%badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=ones(samplesize,tracelength)*NaN;
for i=1:samplesize
    signal=sampletraces(i,stats(i,1):stats(i,2));
    sigstore(i,stats(i,1):stats(i,2))=signal;
    %%% calc degend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sig_height=signal;
    sig_time=(1:length(signal))/length(signal);
    sig_revslope=getslope_reverse(signal,1:slopewindow);
    sig_revslopeshort=getslope_reverse(signal,1:slopewindowshort);
    sig_fwdslope=getslope_forward(signal,1:slopewindow);
    filterdegend=sig_fwdslope-abs(sig_revslope)-2*sig_height+sig_time;
    %filterdegend=sig_fwdslope-abs(sig_revslopeshort)-2*sig_height+sig_time;
    %degend=find(filterdegend==max(filterdegend),1,'first');
    degend=find(filterdegend(slopewindow+1:end)==max(filterdegend(slopewindow+1:end)),1,'first')+slopewindow;
    %degend=find(filterdegend(slopewindowshort+1:end)==max(filterdegend(slopewindowshort+1:end)),1,'first')+slopewindowshort;
%     if degend<=slopewindow
%         continue;
%     elseif sig_fwdslope(degend)<0.1 || sig_height(degend)>0.2
%         badtraces(i)=1;
%         continue;
%     end
    %altstore(i,stats(i,1):stats(i,2))=abs(sig_revslope)+sig_height;
    altstore(i,stats(i,1):stats(i,2))=sig_fwdslope;
    if sig_fwdslope(degend)<0.075 || sig_height(degend)>0.2
        %badtraces(i)=1;
        continue;
    end
    degends(i)=degend;
    if relative==0
        degends(i)=degend+stats(i,1)-1;
    end
    %%% calc degstart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filterdegstart=sig_revslopeshort-abs(sig_fwdslope)-2*sig_height-sig_time;
%     degstart=find(filterdegstart==max(filterdegstart),1,'first');
%     if degstart<=slopewindow
%         continue;
%     elseif sig_revslope(degstart)<0.1 || sig_height(degstart)>0.2
%         badtraces(i)=1;
%         continue;
%     end
    degstart=find(filterdegstart(slopewindowshort+1:end)==max(filterdegstart(slopewindowshort+1:end)),1,'first')+slopewindowshort;
    if sig_revslope(degstart)<0.1 || sig_height(degstart)>0.2
        %badtraces(i)=1;
        continue;
    end
    filterligaseon=-sig_fwdslope+2*sig_height+sig_time;
    ligaseon=find(filterligaseon(1:degstart-1)==max(filterligaseon(1:degstart-1)),1,'last')+1; %estimate the ligase becomes active
    degstarts(i)=degstart;
    if relative==0
        degstarts(i)=degstart+stats(i,1)-1;
    end

    %altstore(i,stats(i,1):stats(i,2))=sig_fwdslope;
end
%badtraces=badtraces>0;

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
for i=1:samplesize
for i=146:168
for i=101:196
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    x=stats(i,1):stats(i,2);
    plot(x,sigstore(i,x));
    hold on;
    %plot(x,mat2gray(altstore(i,x)),'r');
    plot(x,altstore(i,x),'r');
    ylim([-0.1 1]);
    if ~isnan(degstarts(i))
        plot(degstarts(i),sigstore(i,degstarts(i)),'go','markerfacecolor','g','markersize',6);
    end
    if ~isnan(degends(i))
        plot(degends(i),sigstore(i,degends(i)),'ro','markerfacecolor','r','markersize',6);
    end
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