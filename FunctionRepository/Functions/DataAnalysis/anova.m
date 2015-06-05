clear; close all;
%%%%%%%% define conditions %%%%%%%%%%
conditions = {
    'Media',3,1;
    'DMSO',3,2;
    %'Gefitinib',4,1;
    %'U73122',4,3;
    %'Thapsigargin',3,4;
    %'Forskolin',4,5;
    %'Rapamycin',4,6;
    %'U0126',4,7;
    %'LY294002',4,8;
    %'Ruboxistaurin',3,9;
    %'Harmine',4,10;
    %'Nutlin',4,11;
    %'Ponatinib',4,12;
}
data = [];
offset = [];
for i = 1:size(conditions,1)
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    values = angiepostdrugmitosis(row,col,[21 40]);
    data = [data;values(:)];
    offset = [offset;i*ones(length(values),1)];
end
boxplot(axes,data,offset,'labels',conditions(:,1));
set(gca,'Fontsize',14);
h = findobj(gca,'Type','text')
set(h,'FontSize',14)
set(gcf,'color','w');
%title('Drugspike to Mitosis \n(Spike < 4hrs after previous mitosis)');
saveas(gcf,'h:\Downloads\Fig.jpg');
%ttest2(,0.05)
%[p,h] = ranksum(control,drug,'alpha',0.05)