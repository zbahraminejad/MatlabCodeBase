function plottraces_2(data,xvals,xstring,ystring,ylimits,smoothoption,drugspike,tracestats)
figure, hold on;
for i=1:size(data,1)
    plotdata=data(i,:);
    if smoothoption
        realframes=find(~isnan(data(i,:)));
        plotdata(realframes)=smooth(data(i,realframes));
    end
    line(xvals,plotdata,'color','r');
%     line(xvals+0.2,plotdata,'color','b','linestyle','--','linewidth',3);
%     hold on;
%     scatter(tracestats(i,1)/5,plotdata(tracestats(i,1)),100,'markerfacecolor','w','markeredgecolor','b','linewidth',4);
end

line([0 0],ylimits,'color','k','linestyle','--','linewidth',2);

ylim(ylimits);
xlabel(xstring); xlabel('Time (hrs)');
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'D:\Downloads\FigAllLines.jpg');
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