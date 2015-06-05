function plottraceandIF(data,xvals,xstring,ystring,ylimits,smoothoption)
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