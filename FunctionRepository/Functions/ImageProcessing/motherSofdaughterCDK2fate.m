function motherSofdaughterCDK2fate(Slengths,Cdk2inc,Cdk2low,allmotherstats,maxlengthSG2,framesperhr)
%%% compare mother S-phase lengths by Cdk2inc vs low %%%%%%%%%%%%%%%%%%%%%%
Sinc=Slengths(Cdk2inc)/framesperhr; fprintf('Sinc org n=%0.0f\n',numel(Sinc));
Slow=Slengths(Cdk2low)/framesperhr; fprintf('Slow org n=%0.0f\n',numel(Slow));
numnanSinc=sum(allmotherstats(Cdk2inc,5)>=maxlengthSG2 & isnan(Sinc));
numnanSlow=sum(allmotherstats(Cdk2low,5)>=maxlengthSG2 & isnan(Slow));
fprintf('%0.0f traces too long for Sinc measurement\n',numnanSinc);
fprintf('%0.0f traces too long for Slow measurement\n',numnanSlow);
Sinc(allmotherstats(Cdk2inc,5)<maxlengthSG2 | isnan(Sinc))=[];
Slow(allmotherstats(Cdk2low,5)<maxlengthSG2 | isnan(Slow))=[];
Spval=ranksum(Sinc,Slow);
boxplotdata=[Sinc;Slow];
Snuminc=length(Sinc); Snumlow=length(Slow);
meanSinc=round(mean(Sinc)*10)/10;
meanSlow=round(mean(Slow)*10)/10;
stringinc=['CDK2inc: mean=',num2str(meanSinc),'hrs  (n=',num2str(Snuminc),')'];
stringlow=['CDK2low: mean=',num2str(meanSlow),'hrs  (n=',num2str(Snumlow),')'];
offset=[ones(Snuminc,1);2*ones(Snumlow,1)];
figure,boxplot(axes,boxplotdata,offset,'labels',{stringinc,stringlow});
title(['S length comparison: p-val=',num2str(Spval)]);
set(gcf,'color','w','PaperPosition',[0 0 6 4]);
saveas(gcf,'h:\Downloads\Fig_Slengths.jpg');