function [mitosismarks,lastcellmother]=trackmitosis(tracestats,samplecells)
mitosismarks = cell(numel(samplecells),1);
lastcellmother = zeros(numel(samplecells,1));
for i=1:numel(samplecells)
    orgcellid=samplecells(i);
    cellid=orgcellid;
    condition=true;
    mitosisframes = [];
    savecellid = [];
    while condition
        savecellid = [savecellid cellid];
        mitosisframes = [mitosisframes, tracestats(cellid,1)];
        cellid=tracestats(cellid,4);
        condition=~isnan(cellid);
    end
    if numel(savecellid)==1
        lastcellmother(i) = -1;
    else
        lastcellmother(i) = savecellid(2);
    end
    mitosismarks{i} = wrev(mitosisframes);
end