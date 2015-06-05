function scatterandfit(x,y,dotcolor,barcolor,bins)
scatter(x,y,dotcolor,'o');
linespec=['-',barcolor];
binstep=bins(2)-bins(1);
for i=bins
    binsample=y(x>i-binstep/2 & x<=i+binstep/2);
    if isempty(binsample)
        continue;
    end
    [errmean,errstd]=normfit(binsample);
    %errorbar(i,errmean,errstd,color,'-','linewidth',2);
    errorbar(i,errmean,errstd/sqrt(numel(binsample)),linespec,'linewidth',4);
end