function data=incorporatejitters(data,x,y,dims)
fmin=1;fmax=size(data,3);
%%% shift coordinates of each frame to match frame 1 %%%%
for f=fmin:fmax
    data{f}(:,1) = data{f}(:,1)+x(f);
    data{f}(:,2) = data{f}(:,2)+y(f);
end

%%% establish crop borders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=max(x);
xmax=dims(2)+min(x);
ymin=max(y);
ymax=dims(1)+min(y);

%%% eliminate cells falling outside of borders %%%%%%%%%%
for f=fmin:fmax
    outsiders = find(data{f}(:,1)<xmin | data{f}(:,1)>xmax | data{f}(:,2)<ymin | data{f}(:,2)>ymax);
    data{f}(outsiders,:)=[];
end
