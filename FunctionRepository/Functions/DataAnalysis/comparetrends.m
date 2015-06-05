function comparetrends(datatotal,xvals,setb,setr,xstring,ystring,ylimits)
datab=datatotal(setb,:);
datar=datatotal(setr,:);
figure, hold on;
for i=1:size(datab,1)
    line(xvals,datab(i,:),'color','b');
end
for i=1:size(datar,1)
    line(xvals,datar(i,:),'color','r');
end
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\FigCompareLines.jpg');

figure, hold on;
binstep=0.2; %absolute bin width (usu 1frame=0.2hrs)
bincurveshade(xvals,datab,binstep,'b');
bincurveshade(xvals,datar,binstep,'r');
ylim(ylimits);
xlabel(xstring);
ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\FigTrends.jpg');