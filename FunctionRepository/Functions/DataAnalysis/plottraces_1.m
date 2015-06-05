function plottraces_1(data,xvals,xstring,ystring,ylimits,smoothoption)
figure, hold on;
for i=1:size(data,1)
    plotdata=data(i,:);
    if smoothoption
        realframes=find(~isnan(data(i,:)));
        plotdata(realframes)=smooth(data(i,realframes));
    end
    line(xvals,plotdata,'color','b');
end
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
%saveas(gcf,'h:\Downloads\FigAllLines.jpg');

figure, hold on;
binstep=0.2; %absolute bin width (usu 1frame=0.2hrs)
bincurveshade(xvals,data,binstep,'b');
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
%saveas(gcf,'h:\Downloads\FigAllTrends.jpg');
hold off;