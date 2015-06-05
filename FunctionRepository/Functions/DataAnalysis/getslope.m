function theta=getslope(series,step)
%ns=[series,series(end-1:-1:end-step)]; %slope ahead of point
%theta=atan2(ns(1+step:end)-ns(1:end-step),step);
theta=atan2(series(1+step:end)-series(1:end-step),step);
theta=[theta,ones(1,step)*theta(end)];