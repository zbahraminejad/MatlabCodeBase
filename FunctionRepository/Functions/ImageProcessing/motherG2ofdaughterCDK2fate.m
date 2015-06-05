function motherG2ofdaughterCDK2fate(G2lengths,Cdk2inc,Cdk2low,allmotherstats,maxlengthG2,framesperhr)
%%% compare mother G2 lengths by Cdk2inc vs low %%%%%%%%%%%%%%%%%%%%%%%%%%%
G2inc=G2lengths(Cdk2inc)/framesperhr; fprintf('G2inc org n=%0.0f\n',numel(G2inc));
G2low=G2lengths(Cdk2low)/framesperhr; fprintf('G2low org n=%0.0f\n',numel(G2low));
numnanG2inc=sum(allmotherstats(Cdk2inc,5)>=maxlengthG2 & isnan(G2inc));
numnanG2low=sum(allmotherstats(Cdk2low,5)>=maxlengthG2 & isnan(G2low));
fprintf('%0.0f traces too long for G2inc measurement\n',numnanG2inc);
fprintf('%0.0f traces too long for G2low measurement\n',numnanG2low);
G2inc(allmotherstats(Cdk2inc,5)<maxlengthG2 | isnan(G2inc))=[];
G2low(allmotherstats(Cdk2low,5)<maxlengthG2 | isnan(G2low))=[];
G2pval=ranksum(G2inc,G2low);
boxplotdata=[G2inc;G2low];
G2numinc=length(G2inc); G2numlow=length(G2low);
meanG2inc=round(mean(G2inc)*10)/10;
meanG2low=round(mean(G2low)*10)/10;
stringinc=['CDK2inc: mean=',num2str(meanG2inc),'hrs  (n=',num2str(G2numinc),')'];
stringlow=['CDK2low: mean=',num2str(meanG2low),'hrs  (n=',num2str(G2numlow),')'];
offset=[ones(G2numinc,1);2*ones(G2numlow,1)];
figure,boxplot(axes,boxplotdata,offset,'labels',{stringinc,stringlow});
title(['G2 length comparison: p-val=',num2str(G2pval)]);
set(gcf,'color','w','PaperPosition',[0 0 6 4]);
saveas(gcf,'h:\Downloads\Fig_G2lengths.jpg');