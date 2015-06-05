function makeCdk2heatmaps(datasorted,Cdk2start,time,heatmapmax,heatmapmin,markoption,alignoption,xstring,framesperhr)
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
if markoption && alignoption==1
    hold on;
    onlyreal=find(Cdk2start~=0);
    Cdk2start_onlyreal=Cdk2start(onlyreal);
    traceid_onlyreal=traceid(onlyreal);
    scatter((Cdk2start_onlyreal)/framesperhr,traceid_onlyreal,8,'r','*');
elseif markoption && alignoption==2
    hold on;
    scatter(zeros(numtraces,1),1:numtraces,12,'r','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%
xlabel(xstring); ylabel('Individual traces');
set(gcf,'color','w','PaperPosition',[0 0 4 7]);
saveas(gcf,'h:\Downloads\FigCdk2Heatmap.jpg');