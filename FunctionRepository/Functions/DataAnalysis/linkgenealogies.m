function fulltraces=linkgenealogies(tracedata,tracestats)
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% find all non-mothers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motherarray=tracestats(:,2);
mothers=unique(motherarray(~isnan(motherarray)));
nonmothers=find(~ismember(1:totalcells,mothers));
%%% complete all non-mothers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numnonmothers=numel(nonmothers);
fulltraces=ones(numnonmothers,totalframes,totalsignals+2)*NaN;
for i=1:numnonmothers
    cellid=nonmothers(i);
    condition=true;
    while condition
        goodframes=find(~isnan(tracedata(cellid,:,1)));
        fulltraces(i,goodframes,1:totalsignals)=tracedata(cellid,goodframes,:);
        fulltraces(i,goodframes,totalsignals+1)=cellid;  %label with cellid
        fulltraces(i,goodframes,totalsignals+2)=0;       %default no mitosis = 0
        cellid=motherarray(cellid);
        condition=~isnan(cellid);
        if condition
            fulltraces(i,goodframes(1),totalsignals+2)=1;
        end
    end
end