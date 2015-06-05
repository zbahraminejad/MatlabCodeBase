imagepath='G:\Michael\';
experimentpath='20150401-CC-Diff\';
summarydir = [imagepath,experimentpath,'Data Summary\'];
conditionlist = {'Control','DMI','DMI+CDK4-6i'};
numconditions = length(conditionlist);
% 'Control','DMI','Rosi','Control+CDK4-6i','DMI+CDK4-6i','Rosi+CDK4-6i'
% 'Ctrl-NoGem','DMI-NoGem','Rosi-NoGem','Ctrl-NoGem+CDK4-6i','DMI-NoGem+CDK4-6i','Rosi-NoGem+CDK4-6i'
figure(1), hold on
set(gcf,'color','white','units','normalized','outerposition',[0 0 .7 1]);
figure(2), hold on
set(gcf,'color','white','units','normalized','outerposition',[0 0 .7 1]);
figure(3), hold on
set(gcf,'color','white','units','normalized','outerposition',[0 0 .7 1]);
figure(4), hold on
set(gcf,'color','white','units','normalized','outerposition',[0 0 .7 1]);
binrange = 1:.5:10;
binrange2 = 0:5:120;
for cond = 1:numconditions
    load([imagepath,experimentpath,'Data Summary\no CDK2\',conditionlist{cond},'.mat']); mitosistime = allmarkedmitosis; lastmitosis = lastmitosisframe;
    numgated = length(alltraces1(:,1));
    lastvalue = zeros(numgated,1);
    startmedian = zeros(numgated,1);
    endmedian = startmedian;
    for i = 1:numgated
        lastind = find(~isnan(alltraces1(i,:)),1,'last');
        lastvalue(i) = alltraces1(i,lastind);
        tracestart = alltracestats(i,1);
        traceend = alltracestats(i,2);
        startmedian(i) = median(alltraces1(i,tracestart:tracestart+60));
        endmedian(i) = median(alltraces1(i,traceend-60:traceend));
        sortedvalues = sort(alltraces1(i,tracestart:traceend));
    end
    foldchange = endmedian./startmedian;
    sum(foldchange>3.5)/length(foldchange)
          
    mitosistime = [mitosistime{:}]; mitosistime(mitosistime==1) = [];
    mitosistime = mitosistime./5.5;
    lastmitosis = lastmitosis./5.5;
    figure(1)
    histogram(foldchange,binrange,'normalization','pdf');
    figure(2)
    histogram(lastmitosis,binrange2,'normalization','pdf');
    figure(3)
    cdfplot(lastmitosis)
    figure(4)
    cdfplot(foldchange)
end
figure(1)
title({'Fold Change'},'FontName','Arial','FontSize',30);
set(gca,'FontName','Arial','FontSize',25);
xlabel({'Fold change relative to start'},'FontName','Arial','FontSize',25);
% ylabel({'Counts'},'FontName','Arial','FontSize',25);
legend(conditionlist,'FontName','Arial','FontSize',25);
hold off
figure(2)
title({'Last Mitosis Times'},'FontName','Arial','FontSize',30);
set(gca,'FontName','Arial','FontSize',25);
xlabel({'Time (hour)'},'FontName','Arial','FontSize',25);
% ylabel({'Counts'},'FontName','Arial','FontSize',25);
legend(conditionlist,'FontName','Arial','FontSize',25);
hold off
figure(3)
title({'Last Mitosis - CDF'},'FontName','Arial','FontSize',30);
set(gca,'FontName','Arial','FontSize',25);
xlabel({'Mitosis Time (hour)'},'FontName','Arial','FontSize',25);
% ylabel({'Counts'},'FontName','Arial','FontSize',25);
legend(conditionlist,'FontName','Arial','FontSize',25,'Location','SouthEast');
hold off
figure(4)
title({'Fold Change - CDF'},'FontName','Arial','FontSize',30);
set(gca,'FontName','Arial','FontSize',25);
xlabel({'Fold Change'},'FontName','Arial','FontSize',25);
% ylabel({'Counts'},'FontName','Arial','FontSize',25);
legend(conditionlist,'FontName','Arial','FontSize',25,'Location','SouthEast');
hold off




