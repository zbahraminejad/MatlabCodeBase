rows=2:7;
cols=2:6;
sites=1:4;
numrows = length(rows);
numcols = length(cols);
numsites = length(sites);

imagepath='G:\Zahra\';
experimentpath='20150501-Zahra-Pulse\';

datadir=([imagepath,experimentpath,'Data\']);


%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startframe = 1;
endframe = 547;
frames=startframe:endframe;
channelnames={'Cy5_' ,'YFP_' 'YFP_'};

%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
removequiescent = 0;
IFoption = 0;

conditionlist = {'Control','Rosi','DMI','DMI-1x-12hr','DMI-2x-12hr'};
numblocks = 2000;
for col = 1:numcols
%%% Allocate Memory for all trace storage%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alltraces1 = zeros(numblocks,endframe);
alltraces2 = zeros(numblocks,endframe);
alltraces3 = zeros(numblocks,endframe);
alltracestats = zeros(numblocks,4);
allIDs = zeros(numblocks,1);
allmotherstats = zeros(numblocks,1);
allmarkedmitosis = cell(numblocks,1);
allwellID = zeros(numblocks,3);
for row = 1:numrows
for site = 1:numsites

shot=[num2str(rows(row)),'_', num2str(cols(col)), '_', num2str(sites(site))];
wellname = nameandsite(shot);
 rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datafile =[datadir,'tracedata_',shot,'_nolink','.mat'];
% load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
load(datafile,'tracedata','genealogy','jitters');

tic
[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother]=...
    gathertracedata_mz_1(datadir,datafile,shot,motheroption,daughteroption,IFoption);
toc
% badframe = 233;
% [tracedata,jitters]=interpolateframes_MZ(tracedata,jitters,badframe);

[tempcellcount,~]=size(tracestats);
alteredlength  = zeros(tempcellcount,1); alteredstart = alteredlength; gatelength = alteredlength;

for i = 1:tempcellcount
    alteredlength(i) = sum(~isnan(tracedata(i,:,1)));
    alteredstart(i) = find(~isnan(tracedata(i,:,1)),1,'first');
    gatelength(i) = sum(~isnan(tracedata(i,24:end,1)));
end

tracestats(:,1) = alteredstart;
tracestats(:,3) = alteredlength;
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace=380;
badlengths=gatelength<minlengthtrace;

%%% gate degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelPPARg=6; nucareachannel = 3;
noisethresh = 10^6;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=10^10;  %threshold above which max of each trace must be %50
minpos = 0;      %0:anywhere 1:firstframe 2:lastframe
minthresh = 0; %threshold below which min of each trace must be %50

[traces1,badtraces1]=gate_pparg_1_allsites(tracedata,tracestats,...
    noisethresh,channelPPARg,nucareachannel,minthresh,minpos,maxthresh,maxpos);
%%% gate miscellaneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces2=tracedata(:,:,3); %%% Area
traces3=tracedata(:,:,4); % Mass (total h2b intensity median*area);

% orphan = isnan(tracestats(:,4))& ~startexist;
badstart = tracestats(:,1)> 48;
overlap = zeros(tempcellcount,1);
repeatedcelltrace = unique(lastcellmother,'sorted');
repeatedcelltrace(repeatedcelltrace==-1)=[];
overlap(repeatedcelltrace) = 1; overlap(lastcellmother==-1)=0;
%%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces = badlengths | badtraces1 | overlap | badstart;% badtraces2;% | orphan;  
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
traces3=traces3(~badtraces,:);
tracedata=tracedata(~badtraces,:,:);
tracestats=tracestats(~badtraces,:);
markedmitosis = markedmitosis(~badtraces);
IDs=IDs(~badtraces,:);
wellID = zeros(length(IDs),3);
wellID(:,1) = rows(row);
wellID(:,2) = cols(col);
wellID(:,3) = sites(site);
motherstats = motherstats(~badtraces);

%%%%Store trace information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numtraces,~] = size(tracestats);
storeSTARTindex = find(alltracestats(:,1)==0,1,'first');
storeENDindex = storeSTARTindex + (numtraces-1);
lengthpersite = zeros(numrows,1);
lengthpersite(row) = length(traces1(:,1));
% add blocks if not sufficient

if storeENDindex > length(allIDs)
    alltraces1 = [alltraces1; zeros(numblocks,endframe)];
    alltraces2 = [alltraces2; zeros(numblocks,endframe)];
    alltraces3 = [alltraces3; zeros(numblocks,endframe)];
    alltracestats = [alltracestats;zeros(numblocks,4)];
    allIDs = [allIDs;zeros(numblocks,1)];
    allwellID = [allwellID ; zeros(numblocks,3)];
    allmotherstats = [allmotherstats; zeros(numblocks,1)];
    allmarkedmitosis = [allmarkedmitosis;cell(numblocks,1)];
end

% store information 

    alltraces1(storeSTARTindex:storeENDindex,:) = traces1;
    alltraces2(storeSTARTindex:storeENDindex,:) = traces2;
    alltraces3(storeSTARTindex:storeENDindex,:) = traces3;
    alltracestats(storeSTARTindex:storeENDindex,:) = tracestats;
    allIDs(storeSTARTindex:storeENDindex) = IDs;
    allwellID(storeSTARTindex:storeENDindex,:) = wellID;
    allmotherstats(storeSTARTindex:storeENDindex) =  motherstats;
    allmarkedmitosis(storeSTARTindex:storeENDindex) = markedmitosis;
    
end
end

%%%%%%%%%% cut excess blocks from stored matrices and cells %%%%%%%%%%%%%%%
excessblocks = allIDs==0;
alltraces1(excessblocks,:) = [];
alltraces2(excessblocks,:) = [];
alltraces3(excessblocks,:) = [];
alltracestats(excessblocks,:) = [];
allIDs(excessblocks) = [];
allwellID(excessblocks,:) = [];
allmotherstats(excessblocks) = [];
allmarkedmitosis(excessblocks) = [];

numgated=length(alltracestats(:,1));
lastmitosisframe = cellfun(@(x) x(end),allmarkedmitosis);

%%%%%%%%%%%%%%%%%%%%%%%%%save collected data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([imagepath,experimentpath,'Data Summary\',conditionlist{col},'.mat'],'alltraces1','alltraces2','alltraces3',...
    'alltracestats','allIDs','allwellID','allmotherstats','allmarkedmitosis','lastmitosisframe');
end
