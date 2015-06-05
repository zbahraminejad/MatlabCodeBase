function avgtheta=getslope_reverse(series,step)
numsteps=length(step);
thetalength=length(series)-max(step);
theta=zeros(numsteps,thetalength);
for cc=1:numsteps
    i=step(cc);
    temp=atan2(series(1:end-i)-series(1+i:end),i);
    theta(cc,:)=temp(end-thetalength+1:end);
end
avgtheta=mean(theta,1);
avgtheta=[ones(1,max(step))*avgtheta(1),avgtheta];