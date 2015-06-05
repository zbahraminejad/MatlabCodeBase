function plottraces_2(data,xvals,xstring,ystring,xlimits,ylimits,smoothoption,drugspike,tracestats)
figure, hold on;
for i=1:size(data,1)
    plotdata=data(i,:);
    if smoothoption
        realframes=find(~isnan(data(i,:)));
        plotdata(realframes)=smooth(data(i,realframes),5); %default 5
    end
    line(xvals,plotdata,'color','b');
%     line(xvals+0.2,plotdata,'color','b','linestyle','--','linewidth',3);
%     hold on;
%     scatter(tracestats(i,1)/5,plotdata(tracestats(i,1)),100,'markerfacecolor','w','markeredgecolor','b','linewidth',4);
end

line([0 0],ylimits,'color','k','linestyle','--','linewidth',2);

axis([xlimits ylimits]);
set(gca,'fontsize',16);
xlabel(xstring,'fontsize',16); %xlabel('Time (hrs)');
ylabel(ystring,'fontsize',16);
set(gcf,'color','w','PaperPosition',[0 0 4 3]); %default [0 0 4 3] %long [0 0 6 2]
saveas(gcf,'D:\Downloads\FigAllLines.jpg');
saveas(gcf,'D:\Downloads\FigAllLines.eps');
% 
% figure, hold on;
% %hold on;
% binstep=0.2; %absolute bin width (usu 1frame=0.2hrs)
% bincurveshade(xvals,data,binstep,'b');
% ylim(ylimits);
% xlabel(xstring);
% ylabel(ystring);
% %line([0 0],ylimits,'color','r','linewidth',2);
% xlim([-2 4]);
% set(gcf,'color','w','PaperPosition',[0 0 4 4]); %[4 3]
% saveas(gcf,'h:\Downloads\FigAllTrends.jpg');
% hold off;