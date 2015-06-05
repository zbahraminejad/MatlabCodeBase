function [peaklowerbound,peakmid,peakupperbound]=ParameterizeFirstPeak(ImageOld)
TempSeries=ImageOld(:);

%%% get mode of first peak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmax=max(TempSeries);
tmin=min(TempSeries);
nbin=200;tbin=(tmax-tmin)/nbin;
tmin=tmin+tbin/2;tmax=tmax-tbin/2;
[n,xout]=ksdensity(TempSeries,tmin:tbin:tmax);
gp=max([2,ceil(nbin/50)]);      %gp=4
ng=getcurvature(n,gp);          %returns the angle change of intensity histogram over one bin at each point
[~,Ibg00]=min(ng);              %find most convex point of intensity histogram
bg0=xout(Ibg00);

TempSeries2=TempSeries((TempSeries>(bg0-5*gp*tbin))&(TempSeries<(bg0+5*gp*tbin)));  %narrow search for peak
[n_bg,xout_bg]=ksdensity(TempSeries2);
[~,Ibg]=max(n_bg);              %index of peak
peakmid = xout_bg(Ibg);

%%% get upperbound of first peak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
farright = bg0+10*gp*tbin;
[n_th,xout_th]=ksdensity(TempSeries,peakmid:(farright-peakmid)/100:farright);
ng_th=getcurvature(n_th,gp);
pospeak=regionprops(bwlabel(ng_th>0),'PixelIdxList','Centroid','Area');    %same as negpeak, but for positive (concave) curvature
Ith00=zeros(1,size(pospeak,1));SC=Ith00;
for cc=1:size(pospeak,1)
    SC(cc)=sum(ng_th(pospeak(cc).PixelIdxList));                           %sum of curvatures for each pospeak group
    Ith00(cc)=round(pospeak(cc).Centroid(1));                              %index of center of each pospeak group
end
[~,I]=max(SC);
peakupperbound=xout_th(Ith00(I));                                         %I added 10 so the threshold is higher

%%% get lowerbound of first peak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
farleft = bg0-10*gp*tbin;
[n_th,xout_th]=ksdensity(TempSeries,farleft:(peakmid-farleft)/100:peakmid);
ng_th=getcurvature(n_th,gp);
pospeak=regionprops(bwlabel(ng_th>0),'PixelIdxList','Centroid','Area');    %same as negpeak, but for positive (concave) curvature
Ith00=zeros(1,size(pospeak,1));SC=Ith00;
for cc=1:size(pospeak,1)
    SC(cc)=sum(ng_th(pospeak(cc).PixelIdxList));                           %sum of curvatures for each pospeak group
    Ith00(cc)=round(pospeak(cc).Centroid(1));                              %index of center of each pospeak group
end
[~,I]=max(SC);
peaklowerbound=xout_th(Ith00(I));