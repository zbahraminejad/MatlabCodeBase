function [tracedata,daughterstats,motherstats,IFdata]=getdatadaughtermother(datadir,shot,IFoption)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy');
if IFoption
    load([datadir,'IF_',shot,'.mat'],'IFdata');
end
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~IFoption
    samplecells=find(~isnan(genealogy));
    IFdata=[];
else
    samplecells=find(~isnan(genealogy) & ~isnan(IFdata(:,1)));
    IFdata=IFdata(samplecells,:);
end
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=linkancestry(tracedata,tracestats,samplecells);
%%% record mother stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[daughterstats,motherstats,tracedata]=getmotherstats(tracedata,tracestats,samplecells);
end