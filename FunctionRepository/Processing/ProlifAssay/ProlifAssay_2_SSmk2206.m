cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';

load([path,'cellcountmatrix_MK2206_0hr.mat']);
t0=cellcountmatrix;
load([path,'cellcountmatrix_MK2206_24hr.mat']);
t1=cellcountmatrix;
tdiff=t1-t0;
foldchanges=(tdiff./t0)+1;    %normalized difference
mean(foldchanges)

controlfig=figure;
set(controlfig,'color','w','PaperPosition',[0 0 3 3]); %cumulative stairs (3x3)
subaxis(1,1,1,'ML',0.2,'MR',0.2','MT',0.1,'MB',0.1);
notBoxPlot(foldchanges);
set(gca,'xtick',1,'xticklabel','0.1% DMSO');
xlim([0 2]);
ylim([1 2.5]);
ylabel('Fold Change');
saveas(controlfig,'h:\Downloads\FigControl.jpg');
close(controlfig);

cd ..\Processing; %return to this directory