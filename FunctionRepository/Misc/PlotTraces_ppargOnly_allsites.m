% allLength = [];
clear all 
% close all

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
% projectpath='D:\Documents\Projects\';
imagepath='G:\Zahra\';
% imagepath='E:\';
experimentpath='20150501-Zahra-Pulse\';
datadir=([imagepath,experimentpath,'Data Summary\']);


%%%%%%%% Data Selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'DMI','.mat'])
gatetraces = alltraces1./alltraces2;
% alltraces1 = alltraces1./alltraces2;

%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 1;    %0:one trace per plotn 1:all traces in one plot
selectmode = 0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2 = 0;
plotsignal3 = 0;
immunoframe=0;
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
removequiescent = 0;
%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=8;
framesperhr=5.5;
drugspike=17/framesperhr; %114
drugspike2 =93/framesperhr; %386
drugspike3 =155/framesperhr; %386
drugspike4 = 234/framesperhr; %386
startframe = 1;
endframe = 547;%547;
frames=startframe:endframe;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1 = 250000;
ymin2=0; ymax2=3000; %ymax2=1000;

numgated=length(alltracestats(:,1))

lastvalue = zeros(numgated,1);
for i = 1:numgated
    lastind = find(~isnan(alltraces1(i,:)),1,'last');
    lastvalue(i) = alltraces1(i,lastind);
end
% histogram(lastvalue,20)

for i = 1:numgated
   tracestart = alltracestats(i,1);
   traceend = alltracestats(i,2);
   temptrace = gatetraces(i,:);
%    temptrace = alltraces1(i,:);
   startmedian(i) = median(temptrace(tracestart:tracestart+30));
   endmedian(i) = median(temptrace(traceend-30:traceend));
   sumthresh(i) = sum(temptrace(:)>40000);
%    sortedvalues = sort(alltraces2(i,tracestart:traceend));
%     numbmitosis = cellfun(@(x) (length(x)-1),allmarkedmitosis);
end
foldchange = endmedian./startmedian;
highcells = endmedian>200;
disp((sum(highcells)/length(highcells))*100)
% goodlengths = alltracestats(:,3)>450;
gate =  highcells;
% histogram(foldchange,25)
alltraces1 = alltraces1(gate,:) ;
alltraces2 = alltraces2(gate,:) ;
alltraces3 = alltraces3(gate,:) ;
alltracestats = alltracestats(gate,:);
allwellID = allwellID(gate,:);
allIDs = allIDs(gate);
allmarkedmitosis = allmarkedmitosis(gate);
numgated = length(allIDs)
% if numgated>192
%     indexvector = 1:numgated;
%     numgated=192;
%     randomized = sort(datasample(indexvector,numgated,'Replace',false));
%     alltraces1 = alltraces1(randomized,:) ;
%     alltraces2 = alltraces2(randomized,:) ;
%     alltraces3 = alltraces3(randomized,:) ;
%     alltracestats = alltracestats(randomized,:);
%     allwellID = allwellID(randomized,:);
%     allIDs = allIDs(randomized);
%     allmarkedmitosis = allmarkedmitosis(randomized);
% end

%%%data normalization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
% alltraces1=normalizetraces_3(alltraces1,alltracestats,0);
% alltraces2=normalizetraces_3(alltraces2,alltracestats,1);
% alltraces3=normalizetraces_3(alltraces3,alltracestats,1);
%%%%%%%% by median of all max values
% maxval = zeros(numgated,1);
% for ind =  1:numgated
%     maxval(ind) = max(alltraces2(ind,:));
% end
% medianval = median(maxval);
% for ind =  1:numgated
%     alltraces2(ind,:) = alltraces2(ind,:)./medianval;
% end
% for i = 1:numgated
%     temp = allmarkedmitosis{i};
%     temp(temp>endframe) = [];
%     allmarkedmitosis{i} = temp;
% end
% alltraces1 = alltraces1(:,1:endframe);
% select = [60,52];
% numgated =length(select);
% alltraces1 = alltraces1(select,:);
%% 
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=3; %default=2 presentation=4
steps=5; ystep1=round(((ymax1-ymin1)/steps)*10)/10; ystep2=round(((ymax2-ymin2)/steps)*10)/10;
trace2color=[0 0 1]; % blue
% trace2color=[0.8 0.8 0]; % cyan
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
    drugtime2=drugspike2;
    drugplot2=0;
    drugtime3=drugspike3;
    drugplot3=0;
    drugtime4=drugspike4;
    drugplot4=0;
else
    drugtime=0;
    drugplot=drugspike;
    
    drugtime2=0;
    drugplot2=drugspike2;
    
    drugtime3=0;
    drugplot3=drugspike3;
    
    drugtime4=0;
    drugplot4=drugspike4;
end
xtime=frames(1:end)/framesperhr;
% xtime = time/framesperhr;
% xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(IDs,selectin));
end
if ensemble
     figure(2); hold on;
%      subtightplot(numrows,numcols,posmat(row,col),[0.05 0.01],0.02,0.02); hold on;
end
for counter=1:length(selection)
    i=selection(counter);
    if selectmode
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        plotfignum=gcf;
        %set(gca,'Position',[0.03 0.1 0.95 0.8]);
    elseif ~ensemble
        figure(ceil(counter/24));
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(alltraces1(i,:),5); %default 3
    if ~plotsignal2
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth,'Color','b');
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'FontName','Arial','FontSize',16);
    elseif plotsignal2
        ysig2=smoothignorenans(alltraces2(i,:),3);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YColor','k','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color',[0.9 .75 0],'linewidth',tracewidth); % was blue before 'b'
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike>0
        line([drugplot drugplot],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
        line([drugtime drugtime],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
    end
    
    if drugspike2>0
        line([drugplot2 drugplot2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
        line([drugtime2 drugtime2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    
    if drugspike3>0
        line([drugplot3 drugplot3],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
        line([drugtime3 drugtime3],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
    end
    
    if drugspike4>0
        line([drugplot4 drugplot4],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
        line([drugtime4 drugtime4],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    %%%%% mark all mitosis %%%%%%%%%%%%%%%%%%%%
    if length(allmarkedmitosis{i})>1
        plot(xtime(allmarkedmitosis{i}(2:end)),ysig1(allmarkedmitosis{i}(2:end)),'ko','markerfacecolor', 'k','markersize',dotsize);
    else
        plot(xtime(allmarkedmitosis{i}(1)),ysig1(allmarkedmitosis{i}(1)),'ko','markerfacecolor', 'k','markersize',dotsize);
    end
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2); % Ytick spacing default is --> ymin2:ystep2:ymax2
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if plotsignal3
        ysig3=smoothignorenans(alltraces3(i,:),3);
        line(xtime,ysig3,'color','r','DisplayName',num2str(i),'linewidth',tracewidth); % was green before
    end
%     if ensemble
%         title(strcat('row=', num2str(allwellID(i,1)),', col=',num2str(allwellID(i,2)),', site=',num2str(allwellID(i,3))));
%     else
%         title(strcat('Uid: ', num2str(counter),' IDs=',num2str(allIDs(i)),', row=', num2str(allwellID(i,1)),', col=',num2str(allwellID(i,2)),', site=',num2str(allwellID(i,3))));
%     end
%     title(['Unique ID: ', num2str(counter)],'FontName','Arial','FontSize',12)
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        selectmodedisplay(edgemask);
    end
end
fprintf([num2str(numgated),'\n']);



%{
title('Sample Trace: TFEB-EYFP');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 4 3]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}