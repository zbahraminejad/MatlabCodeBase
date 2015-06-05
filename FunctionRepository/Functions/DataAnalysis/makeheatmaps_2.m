function makeheatmaps_2(datasorted,POI,time,heatmapmax,heatmapmin,markoption,xstring,framesperhr,drugspike)
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
if drugspike
    hold on;
    %scatter(ones(numel(POI),1)*drugspike/framesperhr,traceid,4,'k','*');
    scatter(zeros(numel(POI),1),traceid,4,'k','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%
xlabel(xstring); ylabel('Individual traces');
set(gcf,'color','w','PaperPosition',[0 0 4 7]); %[4 7] [4 7]
saveas(gcf,'h:\Downloads\FigHeatmap.jpg');