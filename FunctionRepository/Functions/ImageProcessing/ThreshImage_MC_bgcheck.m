function [th,isbg]=ThreshImage_MC_bgcheck(pixvals,bgthresh,threshmargin)

%%% get background value: search intensity distribution for most convex point %%%%%%%%%%%%%%%%%%%%%%%
tmax=max(pixvals);
tmin=min(pixvals);
nbin=200;tbin=(tmax-tmin)/nbin;
tmin=tmin+tbin/2;tmax=tmax-tbin/2;
[n,xout]=ksdensity(pixvals,tmin:tbin:tmax);
gp=max([2,ceil(nbin/50)]);  %gp=4
ng=getcurvature(n,gp);      %returns the angle change of intensity histogram over one bin at each point
[~,Ibg00]=min(ng);         %find most convex point of intensity histogram
bg0=xout(Ibg00);

plot(xout,n);
plot(xout,ng);

TempSeries2=pixvals((pixvals>(bg0-5*gp*tbin))&(pixvals<(bg0+5*gp*tbin)));  %narrow search for peak
[n_bg,xout_bg]=ksdensity(TempSeries2);
[~,Ibg]=max(n_bg);              %index of peak
bg = xout_bg(Ibg);

if bg<bgthresh
    isbg=true;
else
    isbg=false;
end

%%% get threshold: search upstream of background for most concave point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upperbound = xout(Ibg00)+10*gp*tbin;
lowerbound = bg;
[n_th,xout_th]=ksdensity(pixvals,lowerbound:(upperbound-bg)/100:upperbound);    %shouldn't ng be lowerbound?
ng_th=getcurvature(n_th,gp);
pospeak=regionprops(bwlabel(ng_th>0),'PixelIdxList','Centroid','Area');    %same as negpeak, but for positive (concave) curvature
Ith00=zeros(1,size(pospeak,1));SC=Ith00;
for cc=1:size(pospeak,1)
    SC(cc)=sum(ng_th(pospeak(cc).PixelIdxList));   %sum of curvatures for each pospeak group
    Ith00(cc)=round(pospeak(cc).Centroid(1));       %index of center of each pospeak group
end
[~,I]=max(SC);
th=xout_th(Ith00(I)+threshmargin);  %I added 10 so the threshold is higher
