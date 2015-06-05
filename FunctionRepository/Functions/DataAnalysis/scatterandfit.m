function scatterandfit(x,y,dotcolor,barcolor,stepsize,minbin)
scatter(x,y,dotcolor,'o');
maxbin=ceil(max(x));
if stepsize==2
    modval=mod(minbin,2);
    if mod(maxbin,2)~=modval
        maxbin=maxbin-1;
    end
end
linespec=['-',barcolor];
for i=minbin:stepsize:maxbin
    binsample=y(x>=i-stepsize/2 & x<i+stepsize/2);
    if isempty(binsample)
        continue;
    end
    [mu,std]=normfit(binsample);
    errorbar(i,mu,std,linespec,'linewidth',2);
    %errorbar(i,mu,std/sqrt(numel(binsample)),linespec,'linewidth',4);
end