function avgcurve=getcurvature_array(series,gp)
numgp=length(gp);
for cc=1:numgp
    i=gp(cc);
    series_padded=[series(1+i:-1:1+1),series,series(end-1:-1:end-i)];
    theta1=atan2(series_padded(1+i:end-i)-series_padded(1:end-2*i),i);
    theta2=atan2(series_padded(1+2*i:end)-series_padded(1+i:end-i),i);
    curve_absolute=acos(cos(theta2-theta1));
    curve_signed=series_padded(1:end-2*i)+series_padded(1+2*i:end)-2*series_padded(1+i:end-i);
    curve(cc,:)=curve_absolute.*(curve_signed>0)-curve_absolute.*(curve_signed<0);
end
avgcurve=sum(curve,1);