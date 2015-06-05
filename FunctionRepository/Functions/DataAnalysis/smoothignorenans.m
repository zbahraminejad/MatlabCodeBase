function finaldata=smoothignorenans(data,smoothres)
firstpoint=find(~isnan(data),1,'first');
lastpoint=find(~isnan(data),1,'last');
smoothtemp=smooth(data(firstpoint:lastpoint),smoothres)';
finaldata=data;
finaldata(firstpoint:lastpoint)=smoothtemp;