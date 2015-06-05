function mitosismarks=trackmitosis(tracedata,tracestats,samplecells)
mitosismarks= cell(numel(samplecells),1);
for i=1:numel(samplecells)
    orgcellid=samplecells(i);
    cellid=orgcellid;
    condition=true;
    while condition
        goodframes=find(~isnan(tracedata(cellid,:,1)));
        linkedtracedata(orgcellid,goodframes,:)=tracedata(cellid,goodframes,:);
        cellid=tracestats(cellid,4);
        condition=~isnan(cellid);
    end
end