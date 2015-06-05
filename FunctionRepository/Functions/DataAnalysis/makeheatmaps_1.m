function makeheatmaps_1(datasorted,POI,time,heatmapmax,heatmapmin,markoption,xstring,framesperhr)
numtraces=size(datasorted,1);
datasorted(datasorted>heatmapmax)=heatmapmax;
datasorted(datasorted<heatmapmin)=heatmapmin;
datasorted(isnan(datasorted))=heatmapmin; %make invalid data black
traceid=1:numtraces;
figure;
imagesc(time,traceid,datasorted);
cmap=colormap(jet);
cmap(1,:)=[0 0 0];
colormap(cmap);
%%% overlay markers %%%%%%%%%%%%%%%%%%%%%%%%%
if markoption
    hold on;
    scatter((POI(~isnan(POI)))/framesperhr,traceid(~isnan(POI)),8,'r','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%
xlabel(xstring); ylabel('Individual traces');
set(gcf,'color','w','PaperPosition',[0 0 3 5]); %[4 7]
saveas(gcf,'h:\Downloads\FigHeatmap.jpg');