function drawHoechstEdUgates(Hoechstval,EdUval,xylim,minH,maxH,minE,maxE)
dscatter(Hoechstval,EdUval);
xlabel('Hoechst sumRFU');
ylabel('EdU log2(sumRFU)');
axis(xylim);
rectangle('Position',[minH,minE,maxH-minH,maxE-minE],'EdgeColor','r','linewidth',2);
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
saveas(gcf,'h:\Downloads\Fig.jpg');