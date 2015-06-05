function [ImageMask,th,bg]=ThreshImage_test_higher(ImageOld)

TempSeries=ImageOld(:);

%%% get background value: search intensity distribution for most convex point %%%%%%%%%%%%%%%%%%%%%%%
tmax=max(TempSeries);
tmin=min(TempSeries);
nbin=200;tbin=(tmax-tmin)/nbin;
tmin=tmin+tbin/2;tmax=tmax-tbin/2;
[n,xout]=ksdensity(TempSeries,tmin:tbin:tmax);
gp=max([2,ceil(nbin/50)]);  %gp=4
ng=getcurvature(n,gp);      %returns the angle change of intensity histogram over one bin at each point

[dummy,Ibg00]=min(ng);         %find most convex point of intensity histogram
bg0=xout(Ibg00);

TempSeries2=TempSeries((TempSeries>(bg0-5*gp*tbin))&(TempSeries<(bg0+5*gp*tbin)));  %narrow search for peak
[n_bg,xout_bg]=ksdensity(TempSeries2);
[dummy,Ibg]=max(n_bg);              %index of peak
xfit=xout_bg(Ibg-1:Ibg+1)';yfit=n_bg(Ibg-1:Ibg+1)';
det_n=det([xfit.^2,xfit,[1;1;1]]);
a=det([yfit,xfit,[1;1;1]])/det_n;
b=det([xfit.^2,yfit,[1;1;1]])/det_n;
bg=-b/a/2;                          %refined index of peak --> background intensity

%%% get threshold: search upstream of background for most concave point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upperbound = xout(Ibg00)+10*gp*tbin;
lowerbound = bg;
[n_th,xout_th]=ksdensity(TempSeries,lowerbound:(upperbound-bg)/100:upperbound);    %shouldn't ng be lowerbound?
ng_th=getcurvature(n_th,gp);
pospeak=regionprops(bwlabel(ng_th>0),'PixelIdxList','Centroid','Area');    %same as negpeak, but for positive (concave) curvature
Ith00=zeros(1,size(pospeak,1));SC=Ith00;
for cc=1:size(pospeak,1)
    SC(cc)=sum(ng_th(pospeak(cc).PixelIdxList));   %sum of curvatures for each pospeak group
    Ith00(cc)= max(pospeak(cc).PixelIdxList);       %index of center of each pospeak group
end
[dummy,I]=max(SC);
th=xout_th(Ith00(I));

%%% get mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImageMask=single(ImageOld>th);