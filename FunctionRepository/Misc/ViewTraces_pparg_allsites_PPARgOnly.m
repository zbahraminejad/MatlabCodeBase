%  function ViewTraces
rows=2:7;
cols=2;
sites=1:4;
numrows = length(rows);
numcols = length(cols);
numsites = length(sites);
posmat = 1:(numrows*numcols);
posmat = reshape(posmat,numcols,numrows)';


% allLength = [];

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
% projectpath='D:\Documents\Projects\';
imagepath='G:\Zahra\';
% imagepath='E:\';
experimentpath='20150501-Zahra-Pulse\';

% shot=wellnum2str(row,col,site);
datadir=([imagepath,experimentpath,'Data\']);
% separatedirectories=0;
% if separatedirectories==1
%     rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
%     %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
%     maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
% else
%     rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
%     maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
% %     rawdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
% end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=8;
framesperhr=6;
drugspike=0/framesperhr;
drugspike2 = 0/framesperhr;
startframe = 1;
endframe = 653;
frames=startframe:endframe;
channelnames={'Cy5_' ,'YFP_' 'YFP_'};
edgemask=1; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 0;    %0:one trace per plotn 1:all traces in one plot
selectmode = 0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2 = 0;
plotsignal3 = 0;
IFoption=0; %0:no IFdata 1:IFdata
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1 =250000;
ymin2=0; ymax2=3; %ymax2=1000;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
removequiescent = 0;

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
[tempcellcount,~]=size(tracestats);
alteredlength  = zeros(tempcellcount,1); alteredstart = alteredlength; gatelength = alteredlength;

for i = 1:tempcellcount
    alteredlength(i) = sum(~isnan(tracedata(i,:,1)));
    alteredstart(i) = find(~isnan(tracedata(i,:,1)),1,'first');
    gatelength(i) = sum(~isnan(tracedata(i,24:end,1)));
end

tracestats(:,1) = alteredstart;
tracestats(:,3) = alteredlength;
% histogram(tracestats(:,3),100)
% hist(tracedata(:,end,7),50)
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
% figure
% reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace=400;
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

%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
% traces1=normalizetraces_3(traces1,tracestats,1);
% traces2=normalizetraces_3(traces2,tracestats,1);
% traces3=normalizetraces_3(traces3,tracestats,0);

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
save([imagepath,experimentpath,'Data\Corrected\',conditionlist{cols(col)},'.mat'],'alltraces1','alltraces2','alltraces3',...
    'alltraces4','alltraces5','alltracestats','allIDs','allwellID',...
    'allmotherstats','allmarkedmitosis','lastmitosisframe');
end
%%
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
steps=5; ystep1=round(((ymax1-ymin1)/steps)*10)/10; ystep2=round(((ymax2-ymin2)/steps)*10)/10;
trace2color=[0 0 1]; % blue
% trace2color=[0 1 1]; % cyan
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
    drugtime2=drugspike2;
    drugplot2=0;
else
    drugtime=0;
    drugplot=drugspike;
    drugtime2=0;
    drugplot2=drugspike2;
end
xtime=frames/framesperhr;
% xtime = time/framesperhr;
% xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(IDs,selectin));
end
if ensemble
     figure; hold on;
%      subtightplot(numrows,numcols,posmat(row,col),[0.05 0.01],0.02,0.02); hold on;
end
for counter=1:length(selection)
    i=selection(counter);
    if selectmode
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        plotfignum=gcf;
        %set(gca,'Position',[0.03 0.1 0.95 0.8]);
    elseif ~ensemble
        figure(ceil(counter/24));
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(alltraces1(i,:),12); %default 3
    if ~plotsignal2
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
    elseif plotsignal2
        ysig2=smoothignorenans(alltraces2(i,:),6);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YColor','k','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color',[0.9 .75 0],'linewidth',tracewidth); % was blue before 'b'
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike>0
        line([drugplot drugplot],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
        line([drugtime drugtime],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    
    if drugspike2>0
        line([drugplot2 drugplot2],[ymin1 ymax1],'Color','k','linewidth',1,'linestyle','-');
        line([drugtime2 drugtime2],[ymin1 ymax1],'Color','k','linewidth',1,'linestyle','-');
    end
    %%%%% mark all mitosis %%%%%%%%%%%%%%%%%%%%
    if ~isnan(alltracestats(i,4))
        plot(xtime(allmarkedmitosis{i}(2:end)),ysig1(allmarkedmitosis{i}(2:end)),'ro','markerfacecolor', 'g','markersize',dotsize);
    end
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:1:ymax2); % Ytick spacing default is --> ymin2:ystep2:ymax2
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if plotsignal3
        ysig3=smoothignorenans(alltraces3(i,:),6);
        line(xtime,ysig3,'color','r','DisplayName',num2str(i),'linewidth',tracewidth); % was green before
    end
    if ensemble
        title(strcat('row=', num2str(allwellID(i,1)),', col=',num2str(allwellID(i,2)),', site=',num2str(allwellID(i,3))));
    else
        title(strcat('IDs=',num2str(allIDs(i)),', row=', num2str(allwellID(i,1)),', col=',num2str(allwellID(i,2)),', site=',num2str(allwellID(i,3))));
    end
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        selectmodedisplay(edgemask);
    end
end
fprintf([num2str(numgated),'\n']);



%{
title('Sample Trace: TFEB-EYFP');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 4 3]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}