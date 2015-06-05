
% allLength = [];

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
% projectpath='D:\Documents\Projects\';
imagepath='G:\Michael\';
% imagepath='E:\';
experimentpath='20150401-CC-Diff\';
datadir=([imagepath,experimentpath,'Data Summary\Corrected\']);


%%%%%%%% Data Selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'Rosi','.mat'])
% temphold = alltraces1;
alltraces2 = alltraces3;
% alltraces3 = temphold;

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
drugspike=114/framesperhr; %114
drugspike2 = 386/framesperhr; %386
startframe = 1;
endframe = 653;
frames=startframe:endframe;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=10^3; ymax1 = 2.5*10^5;
ymin2=0; ymax2=3000; %ymax2=1000;

numgated=length(alltracestats(:,1))
% firstmitosis = 
% for tracenum = 1:numgated
%    firstmitosis 
% end
% [~,sortedind] = sort(lastmitosisframe,'ascend');
% sorttraces = alltraces1;
% sorttraces = alltraces1(sortedind,:);
% figure,imagesc(sorttraces,[0 25000]);
% firstmitosis = cellfun(@(x) x(2), allmarkedmitosis);
% [time,alltraces1] = aligntraces_5(alltraces1,firstmitosis,alltracestats,allmotherstats,daughteroption);
% [~,alltraces2] = aligntraces_5(alltraces2,firstmitosis,alltracestats,allmotherstats,daughteroption);

% 
% [~,sortedlastindex] = sort(lastmitosisframe);
% [~,sortedfirstindex] = sort(firstmitosis);
% 
% figure,imagesc(alltraces1(sortedlastindex,:),[0 100000]); title({'Sorted by Last Mitosis'},'FontName','Arial','FontSize',25);
% figure,imagesc(alltraces1(sortedfirstindex,:),[0 100000]); title({'Sorted by First Mitosis'},'FontName','Arial','FontSize',25);

lastvalue = zeros(numgated,1);
for i = 1:numgated
    lastind = find(~isnan(alltraces1(i,:)),1,'last');
    lastvalue(i) = alltraces1(i,lastind);
end
% histogram(lastvalue,20)

lastmitosis = cellfun(@(x) x(end), allmarkedmitosis);
g2exit = zeros(numgated,1);
for cell = 1:numgated
   temptrace = smoothignorenans(alltraces3(cell,lastmitosis(cell):end),3);
   sumthresh = sum(temptrace>200);
   if sumthresh > 6
   g2exit(cell) = 1;
   end
end
disp(sum(g2exit)/length(g2exit)*100)
g2exit = logical(g2exit);

startmedian = zeros(numgated,1);
sumthresh = zeros(numgated,1);
endmedian = startmedian;
% numbmitosis = zeros(numgated,1);
for i = 1:numgated
   tracestart = alltracestats(i,1);
   traceend = alltracestats(i,2);
   startmedian(i) = median(alltraces1(i,tracestart:tracestart+33));
   endmedian(i) = median(alltraces1(i,traceend-33:traceend));
   sumthresh(i) = sum(alltraces1(i,:)>40000);
%    sortedvalues = sort(alltraces2(i,tracestart:traceend));
%     numbmitosis = cellfun(@(x) (length(x)-1),allmarkedmitosis);
end
foldchange = endmedian./startmedian;
highcells = foldchange>4 & sumthresh>11;
disp((sum(highcells)/length(highcells))*100)
% goodlengths = alltracestats(:,3)>450;
gate =  highcells;
% histogram(foldchange,25)
alltraces1 = alltraces1(gate,:) ;
alltraces2 = alltraces2(gate,:) ;
alltraces3 = alltraces3(gate,:) ;
alltraces4 = alltraces4(gate,:) ;
alltraces5 = alltraces5(gate,:) ;
alltracestats = alltracestats(gate,:);
allwellID = allwellID(gate,:);
allIDs = allIDs(gate);
allmarkedmitosis = allmarkedmitosis(gate);
numgated = length(allIDs)
% if numgated>150
%     indexvector = 1:numgated;
%     numgated=150;
%     randomized = sort(datasample(indexvector,numgated,'Replace',false));
%     alltraces1 = alltraces1(randomized,:) ;
%     alltraces2 = alltraces2(randomized,:) ;
%     alltraces3 = alltraces3(randomized,:) ;
%     alltracestats = alltracestats(randomized,:);
%     allwellID = allwellID(randomized,:);
%     allIDs = allIDs(randomized);
%     allmarkedmitosis = allmarkedmitosis(randomized);
% end
% for cell =1:numgated
%    alltraces1(cell,:) = smoothignorenans(alltraces1(cell,:),6); 
% end    
% alltraces1 = diff(alltraces1,1,2);

% histogram(lastmitosisframe(highcells),'normalization','pdf')
% hold on
% histogram(lastmitosisframe(~highcells),'normalization','pdf')
% hold off
% cdfplot(lastmitosisframe(highcells))
% hold on
% cdfplot(lastmitosisframe(~highcells))
% hold off
% set(gca,'FontName','Arial','FontSize',36)
% xlabel({'Frame Number'},'FontName','Arial','FontSize',36)
% title({'Last Mitosis Frame - CDF'},'FontName','Arial','FontSize',36)
% legend({'PPARg high','PPARg low'},'FontName','Arial','FontSize',36)

%%%data normalization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
% alltraces1=normalizetraces_3(alltraces1,alltracestats,1);
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

%%%%% to look at slope differences in geminin traces%%%%%%%%%%%%%%%%%%%
% for cell = 1:numgated
%     alltraces2(cell,:) = smoothignorenans(alltraces1(cell,:),3);
% end
% alltraces1 = diff(alltraces2,1,2);
% alltraces2 = alltraces2(:,2:end);
% alltraces3 = alltraces3(:,2:end);
%% 
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=1; %default=2 presentation=4
steps=5; ystep1=round(((ymax1-ymin1)/steps)*10)/10; ystep2=round(((ymax2-ymin2)/steps)*10)/10;
trace2color=[0 0 1]; % blue
% trace2color=[0 1 1]; % cyan
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
    drugtime2=drugspike2;
    drugplot2=0;
else
    drugtime=0;
    drugplot=drugspike;
    drugtime2=0;
    drugplot2=drugspike2;
end
xtime=frames(1:end)/framesperhr;
% xtime = time/framesperhr;
% xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(IDs,selectin));
end
if ensemble
     figure(1); hold on;
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
    ysig1=smoothignorenans(alltraces1(i,:),6); %default 3
    if ~plotsignal2
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth,'Color','r');
        axis([xtime(1) xtime(end) ymin1 ymax1]);
%         set(gca,'FontName','Arial','FontSize',20);
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
        line([drugtime drugtime],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    
    if drugspike2>0
        line([drugplot2 drugplot2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
        line([drugtime2 drugtime2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
    end
    %%%%% mark all mitosis %%%%%%%%%%%%%%%%%%%%
%     if ~isnan(alltracestats(i,4))
        plot(xtime(allmarkedmitosis{i}(1:end)),ysig1(allmarkedmitosis{i}(1:end)),'ko','markerfacecolor', 'k','markersize',dotsize);
%     end
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
    if ensemble
        title(strcat('row=', num2str(allwellID(i,1)),', col=',num2str(allwellID(i,2)),', site=',num2str(allwellID(i,3))));
    else
        title(strcat('IDs=',num2str(allIDs(i)),', row=', num2str(allwellID(i,1)),', col=',num2str(allwellID(i,2)),', site=',num2str(allwellID(i,3))));
    end
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