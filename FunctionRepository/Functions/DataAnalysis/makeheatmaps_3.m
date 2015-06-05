function makeheatmaps_3(datasorted,firstmarker,secondmarker,time,heatmapmax,heatmapmin,xstring,framesperhr)
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
if ~isempty(firstmarker)
    hold on;
    scatter((firstmarker(~isnan(firstmarker)))/framesperhr,traceid(~isnan(firstmarker)),8,'r','*');
end
if ~isempty(secondmarker)
    hold on;
    scatter((secondmarker(~isnan(secondmarker)))/framesperhr,traceid(~isnan(secondmarker)),8,'g','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%
xlabel(xstring); ylabel('Individual traces');
set(gcf,'color','w','PaperPosition',[0 0 4 7]); %[4 7] [4 7]
saveas(gcf,'h:\Downloads\FigHeatmap.jpg');