function plottraces(data,xvals,xstring,ystring,ylimits)
figure, hold on;
for i=1:size(data,1)
    line(xvals,data(i,:),'color','b');
end
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\FigAllLines.jpg');