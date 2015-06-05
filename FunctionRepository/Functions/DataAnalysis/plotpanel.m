function plotpanel(data,xvals,ylimits,smoothoption)
for i=1:size(data,1)
    plotdata=data(i,:);
    if smoothoption
        realframes=find(~isnan(data(i,:)));
        plotdata(realframes)=smooth(data(i,realframes));
    end
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    plot(xvals,plotdata,'color','b');
    ylim(ylimits);
end