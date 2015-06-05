cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Immunofluorescence\2012-12-11 PD0325901 Dose Response\';
experiment = 'pFra-1-488 pErk-594\'; %row 2,4,6  %note: row 6 is no good
%experiment = 'pErk-488 pAkt-594\';   %row 3,5,7
rows = [2 4];

doses = {'0nM','0.01nM','0.03nM','0.1nM','0.3nM','1nM','3nM','10nM','30nM','100nM';};

arraycellcount=zeros(length(rows),10);
arrayred=arraycellcount;
arrayfitc=arraycellcount;
for i=1:length(rows)
    load([path,experiment,'R',num2str(rows(i)),'_alldata_nodapi.mat'],'totalcellcount','totalmedianred','totalmedianfitc');
    arraycellcount(i,:)=mean(totalcellcount(1:10,:),2);
    arrayred(i,:)=mean(totalmedianred(1:10,:),2)';
    arrayfitc(i,:)=mean(totalmedianfitc(1:10,:),2)';
end
%meancc=mean(arraycellcount);
meanred=mean(arrayred);
stdred=std(arrayred,0,2);
meanfitc=mean(arrayfitc);
stdfitc=std(arrayfitc,0,2);

doseresponsefig=figure;
set(doseresponsefig,'color','w','PaperPosition',[0 0 18 12]); %cumulative stairs (3x3)

numdoses=length(doses);
subaxis(2,1,1,'ML',0.07,'MR',0.02,'MT',0.05,'MB',0.1,'SH',0.04,'SV',0.07);
plot(1:numdoses,meanred,'-mo','linewidth',2,'markeredgecolor','k','markerfacecolor',[.49 1 .63],'markersize',10);
set(gca,'xtick',1:numdoses,'xticklabel',doses);
ylim([2 4.5]);
title('pErk','fontsize',16);
subaxis(2,1,2,'ML',0.07,'MR',0.02,'MT',0.05,'MB',0.1,'SH',0.04,'SV',0.07);
plot(1:numdoses,meanfitc,'-mo','linewidth',2,'markeredgecolor','k','markerfacecolor',[.49 1 .63],'markersize',10);
set(gca,'xtick',1:numdoses,'xticklabel',doses);
ylim([2 4.5]);
title('pErk','fontsize',16);

saveas(doseresponsefig,'h:\Downloads\Fig.jpg');
close(doseresponsefig);
cd ..\Processing; %return to this directory