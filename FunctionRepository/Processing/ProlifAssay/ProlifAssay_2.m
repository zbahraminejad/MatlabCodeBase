cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\2013-04-18_DoseResponsePanel\';

drugs = {
    'PD032591';
    'Gefitinib';
    'GDC-0941';
    'Sorafenib';
    'Ponatinib';
    'U73122';
    'Torin';
    'Nutlin';
    'Imatinib';
    'MK2206';
    'PD0332991';
    'DAGKI-2';
    };

concentrations = {
    'N/A','0.1nM','0.3nM','1nM','3nM','10nM','30nM','100nM';
    '3nM','10nM','30nM','100nM','300nM','1uM','3uM','10uM';
    'N/A','0.1nM','0.3nM','1nM','3nM','10nM','30nM','100nM';
    '0.3nM','1nM','3nM','10nM','30nM','100nM','300nM','1uM'
    'N/A','0.1nM','0.3nM','1nM','3nM','10nM','30nM','100nM';
    '3nM','10nM','30nM','100nM','300nM','1uM','3uM','10uM';
    '0.03nM','0.1nM','0.3nM','1nM','3nM','10nM','30nM','100nM';
    '3nM','10nM','30nM','100nM','300nM','1uM','3uM','10uM';
    '3nM','10nM','30nM','100nM','300nM','1uM','3uM','10uM';
    '3nM','10nM','30nM','100nM','300nM','1uM','3uM','10uM';
    '0.3nM','1nM','3nM','10nM','30nM','100nM','300nM','1uM'
    '0.3nM','1nM','3nM','10nM','30nM','100nM','300nM','1uM'
    };

load([path,'cellcountmatrix_0hr.mat']);
t0=sum(cellcountmatrix,3);
load([path,'cellcountmatrix_24hr.mat']);
t1=sum(cellcountmatrix,3);
tdiff=t1-t0;
tdiffnorm=(tdiff./t0)+1;    %normalized difference

controlfig=figure;
set(controlfig,'color','w','PaperPosition',[0 0 3 3]); %cumulative stairs (3x3)
subaxis(1,1,1,'ML',0.2,'MR',0.2','MT',0.1,'MB',0.1);
controlvals=tdiffnorm(1,[1 3 5]);
notBoxPlot(controlvals);
set(gca,'xtick',1,'xticklabel','0.1% DMSO');
xlim([0 2]);
%ylim([1 5]);
ylim([1 2.5]);
ylabel('Fold Change');
saveas(controlfig,'h:\Downloads\FigControl.jpg');
close(controlfig);
linemean=mean(controlvals);
linestdhigh=linemean+std(controlvals);
linestdlow=linemean-std(controlvals);

doseresponsefig=figure;
set(doseresponsefig,'color','w','PaperPosition',[0 0 18 12]); %cumulative stairs (3x3)
for colidx=1:12
    drug=char(drugs(colidx));
    values=tdiffnorm(:,colidx)';
    doses=concentrations(colidx,:);
    if strcmp(char(doses(1)),'N/A')
        doses=doses(2:end);
        values=values(2:end);
    end
    numdoses=length(doses);
    subaxis(3,4,colidx,'ML',0.07,'MR',0.02,'MT',0.05,'MB',0.1,'SH',0.04,'SV',0.07);
    plot(1:numdoses,values,'-mo','linewidth',2,'markeredgecolor','k','markerfacecolor',[.49 1 .63],'markersize',10);
    set(gca,'xtick',1:numdoses,'xticklabel',doses);
    %ylim([1 5]);
    ylim([1 2.5]);
    hold on; 
    line(1:numdoses,linemean*ones(1,numdoses),'linestyle','--','color','r');
    line(1:numdoses,linestdhigh*ones(1,numdoses),'linestyle','--','color','r');
    line(1:numdoses,linestdlow*ones(1,numdoses),'linestyle','--','color','r');
    title(drug,'fontsize',16);
end
saveas(doseresponsefig,'h:\Downloads\Fig.jpg');
close(doseresponsefig);
cd ..\Processing; %return to this directory