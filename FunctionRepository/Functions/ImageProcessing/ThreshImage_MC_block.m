function bg=ThreshImage_MC_block(ImageOld)
TempSeries=ImageOld(:);

%%% get background value: search intensity distribution for most convex point %%%%%%%%%%%%%%%%%%%%%%%
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
bg = xout_bg(Ibg);
end