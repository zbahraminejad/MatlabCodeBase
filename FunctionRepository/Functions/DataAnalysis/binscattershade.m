function bincurveshade(x,y,stepsize,minbin,maxbin,colorchar)
if stepsize==2
    modval=mod(minbin,2);
    if mod(maxbin,2)~=modval
        maxbin=maxbin-1;
    end
end
bins=minbin:stepsize:maxbin;
numbins=length(bins);
mu=ones(numbins,1)*NaN; std=ones(numbins,1)*NaN;
for i=1:numbins
    b=bins(i);
    binsample=y(x>=b-stepsize/2 & x<b+stepsize/2);
    if isempty(binsample)
        continue;
    end
    %[mu(i),std(i)]=normfit(binsample);
    mu(i)=nanmean(binsample);
    std(i)=nanstd(binsample);
end
binfill=[bins fliplr(bins)];
if strcmp(colorchar,'r')
    colorcode=[1 0 0];
elseif strcmp(colorchar,'g')
    colorcode=[0 1 0];
elseif strcmp(colorchar,'b')
    colorcode=[0 0 1];
end
fadedcolor=colorcode+.7*(1-colorcode);
line(bins',mu,'color',colorcode,'linewidth',2);
fill(binfill',[mu'+std' fliplr(mu'-std')],fadedcolor,'edgecolor',fadedcolor,'FaceAlpha', 0.4);