clear all 
% close all
% rowmat = [2,3,4;5,6,7];
rows= 2:7; cols=[4,7];
% rows = 2:4;
sites = 1:2;
numrows = length(rows);
numsites = length(sites);
numcols = length(cols);
global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
imagepath = 'E:\Michael\';
experimentpath='20150320-Diff Stimulus on CC\';

datadir=([imagepath,experimentpath,'Data\']);
% separatedirectories=0;
% if separatedirectories==1
% rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
% %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
% maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
% else
% %     rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
% maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
% rawdir = [imagepath,experimentpath,'Raw\',wellname,shot,'_'];
% end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
framesperhr=6;
drugspike=100/framesperhr;
startframe = 1; endframe = 199;
frames=startframe:endframe; %319
channelnames={'CFP_' 'YFP_' 'mCherry_'};
edgemask=0; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble=0;    %0:one trace per plotn 1:all traces in one plot
selectmode= 0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2=0;
plotsignal3=0;
IFoption=0; %0:no IFdata 1:IFdata
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1 = 3;
ymin2=0; ymax2=1.5; %ymax2=1000;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=2; %0:no gating 1:mothers (traces that end in mitosis) 2:no mothers (traces that don't end in mitosis)
daughteroption=1; %0:no gating 1:daughters (traces that start with mitosis) 2:no daughters (traces that don't start with mitosis)
quiescentanalysis=0;
removequiescent = 0; % 0 for no gating, 1 to look for cells that stay quiescent the whole time;
%%% Allocate Memory for all trace storage%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numblocks = 5000;
alltraces1 = zeros(numblocks,endframe);
alltraces2 = zeros(numblocks,endframe);
alltraces3 = zeros(numblocks,endframe);
alltracestats = zeros(numblocks,4);
allIDs = zeros(numblocks,1);
allmotherstats = zeros(numblocks,1);
for col = 1:numcols
for row = 1:numrows
    for site = 1:numsites
    shot=[num2str(rows(row)),'_', num2str(cols(col)), '_', num2str(sites(site))];
    wellname = nameandsite(shot);
    rawdir = [imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    %%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    datafile = [datadir,'tracedata_',shot,'_nolink','.mat'];
    load(datafile,'tracedata','genealogy','jitters');
    [tracedata,tracestats,motherstats,IFdata,IDs]=gathertracedata_3(datadir,datafile,shot,motheroption,daughteroption,IFoption);
    %%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
    % reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
    %%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minlengthtrace=24; % THIS IS NOT TOTAL LENGTH SINCE ANCESTRY OF THE CELL IS LINKED AND TOTAL LENGTH IS LONGER FOR MITOTIC CELLS
    badlengths=tracestats(:,3)<minlengthtrace;
    %%% gate degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    channelGem=10;
    maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
    maxthresh=40;  %threshold above which max of each trace must be %50
    minpos=0;      %0:anywhere 1:mothertrace 2:daughtertrace
    minthresh=20; %threshold below which min of each trace must be %50
    % [traces1,badtraces1]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
    [traces2,badtraces2]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
    %%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucchannel=6; cytochannel=7;
    nucthreshoption=0;
    %0:normalize each trace by percentile of entire trace
    %1:normalize by first frame of trace (regardless of whether mitosis or not)
    nucthresh=20;    %threshold for nuc intensity according to nucthreshoption

    motherthresh=0;   %threshold for max DHB ratio of mother. default gating=1.0.  0:ignore
    noisethresh=0.4;  %threshold for max positive jump in cyto/nuc ratio
    [traces1,badtraces1]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis,removequiescent);
    % [traces2,badtraces2]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis,onlyquiescent);
    %%% gate by IF data %%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     staindata = IFdata(:,8);
%     thresholdsig = 0;
%     badIF = (staindata<thresholdsig);
%     orphan = isnan(tracestats(:,4)) & tracestats(:,3)<120;
    
%%%%%% gate miscellaneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    traces3=tracedata(:,:,4); % Area is 3, Mass(total h2b intensity) is 4
    % traces2 = tracedata(:,:,7);
    % traces3 = tracedata(:,:,8);
%     badstart = traces1(:,1)>0.4;
threshgate = zeros(length(tracestats(:,1)),1);
   for i = 1:length(tracestats(:,1))
        threshgate(i) = traces1(i,(tracestats(i,1)))>0.5;
   end
    %%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp(length(tracestats)-sum(badlengths))
    badtraces= badlengths | badtraces1;% |threshgate;  %| orphan; %badtraces1 ;% | badtraces2;
    traces1=traces1(~badtraces,:);
    traces2=traces2(~badtraces,:);
    traces3=traces3(~badtraces,:);
    tracedata=tracedata(~badtraces,:,:);
    tracestats=tracestats(~badtraces,:);
    IDs=IDs(~badtraces,:);
    motherstats = motherstats(~badtraces);
%     staindata = staindata(~badtraces);

    %%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %0:normalize each trace by max of its own trace
    %1:normalize by median of max of all traces if min>max*0.5
    %2:normalize each trace by value at first frame after mitosis
    % traces1=normalizetraces_3(traces1,tracestats,1);
    traces2=normalizetraces_3(traces2,tracestats,0);
    % traces3=normalizetraces_3(traces3,tracestats,0);

%%%%Store trace information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numtraces,~] = size(tracestats);
storeSTARTindex = find(alltracestats(:,1)==0,1,'first');
storeENDindex = storeSTARTindex + (numtraces-1);

% add blocks if not sufficient
if storeENDindex > length(allIDs)
    alltraces1 = [alltraces1; zeros(numblocks,endframe)];
    alltraces2 = [alltraces2; zeros(numblocks,endframe)];
    alltraces3 = [alltraces3; zeros(numblocks,endframe)];
    alltracestats = [alltracestats;zeros(numblocks,4)];
    allIDs = [allIDs;zeros(numblocks,1)];
    allmotherstats = [allmotherstats; zeros(numblocks,1)];
end

% store information 

    alltraces1(storeSTARTindex:storeENDindex,:) = traces1;
    alltraces2(storeSTARTindex:storeENDindex,:) = traces2;
    alltraces3(storeSTARTindex:storeENDindex,:) = traces3;
    alltracestats(storeSTARTindex:storeENDindex,:) = tracestats;
    allIDs(storeSTARTindex:storeENDindex) = IDs;
    allmotherstats(storeSTARTindex:storeENDindex) =  motherstats;
    end
 
end
end
%%%%%%%%%% cut excess blocks from stored matrices and cells %%%%%%%%%%%%%%%
excessblocks = allIDs==0;
alltraces1(excessblocks,:) = [];
alltraces2(excessblocks,:) = [];
alltraces3(excessblocks,:) = [];
alltracestats(excessblocks,:) = [];
allIDs(excessblocks) = [];
allmotherstats(excessblocks) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(7),hist(allstaindata,20)
% title('stain data histogram')
% maxvalue = max(alltraces1,[],2);
% figure(8), hist(maxvalue,30);
% title('hist of max of all traces');
% sum(maxvalue>1.2)/length(maxvalue)

[numcells,~] = size(alltraces1);
% endvalue = zeros(numcells,1);
quiescent = zeros(numcells,1);
for i = 1:numcells
    cycling = sum(alltraces1(i,alltracestats(i,1):end)>.7);
    maxval = max(alltraces1(i,alltracestats(i,1)+12:end));
%     cycling = sum(alltraces1(i,1:end)>.7);
%     maxval = max(alltraces1(i,1:end));
    if cycling >10 || maxval > 1.2
        quiescent(i) = 0;
    else
        quiescent(i) = 1;
    end
end
quiescent = logical(quiescent);
sum(quiescent)/length(quiescent)
[sortedq,index] = sort(quiescent);
cdk2traces = alltraces1;
tracestats = alltracestats;

% [time,alltraces1] = aligntraces_5(alltraces1,alltracestats(:,1),alltracestats,allmotherstats,daughteroption);


[~,sortedmitosisindices]  = sort(alltracestats(:,1),'ascend');
cdk2traces = cdk2traces(sortedmitosisindices,:);
figure(1),imagesc(cdk2traces,[0 2]); colorbar;
% % [~,ind]  = sort(allstaindata,'ascend');
% ind = enterindex;
% % ind = [qind; cind];
% maptrace = maptrace(ind,:);
% maptracestats = alltracestats(cind,:);
% figure(20),imagesc(maptrace, [0 3])
% xlim([temptime(1) temptime(end)])

% [numcells,~] = size(alltraces1);
% endvalue = zeros(numcells,1);
% startind = find(time==0);
% for i = 1:numcells
%     tempsignal = alltraces1(i,startind:end);
%     templength = sum(tempsignal>1);
%     if templength>12
%         endvalue(i)=1;
%     end
% end
% endvalue = logical(endvalue);
% cycling = sum(endvalue);
% alltraces1(~endvalue,:)=[];
% alltracestats(~endvalue,:)=[];
% cycling
% numcells
%% 
% imagesc(alltraces1, [0 3])
numgated=size(alltraces1,1);
if numgated>192
    numgated=120;
end

set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=3; %default=2 presentation=4
steps=5; ystep1=round(((ymax1-ymin1)/steps)*10)/10; ystep2=round(((ymax2-ymin2)/steps)*10)/10;
trace2color=[1 0 0];
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
else
    drugtime=0;
    drugplot=drugspike;
end
xtime = frames/framesperhr;
% xtime=time/framesperhr;
% xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(allIDs,selectin));
end
if ensemble
    figure(3); hold on;
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
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
%         figure(ceil(counter/6));
%         subaxis(2,3,mod(counter-1,6)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.04); %5x4  ML, MR, MT, MB control outter borders and SH controls spacing between subplots
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(alltraces1(i,:),12); %default 3
    if ~plotsignal2
        if quiescent(counter) == 1
            linecolor = [1 0 0];
        elseif quiescent(counter) <0
            linecolor = [0 1 0];
        else
            linecolor = [0 0 1];
        end
%         linecolor = [0 0 1];
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth,'Color',linecolor); % original
%         line(xtime,ysig1,'DisplayName','Control','linewidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        if ensemble
            set(gca,'FontSize',36,'FontName','Arial');
        end
    elseif plotsignal2
        ysig2=smoothignorenans(alltraces2(i,:),6);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike>0
        line([drugplot drugplot],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
        %line([drugtime drugtime],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    %%%%% mark mitosis (only the current cell's birth) %%%%%%%%%%%%%%%%%%%%
%     if ~isnan(alltracestats(i,4))
%         plot(xtime(alltracestats(i,1)),ysig1(alltracestats(i,1)),'ro','markerfacecolor', 'g','markersize',dotsize);
%     end
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2);
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if plotsignal3
        ysig3=smoothignorenans(traces3(i,:),3);
        line(xtime,ysig3,'color','g','DisplayName',num2str(i),'linewidth',tracewidth);
    end
%     title(num2str(allIDs(i)));
%     title('Rosi 1 uM');
%     %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if selectmode
%         selectmodedisplay(edgemask);
    end

fprintf([num2str(numgated),'\n']);
%{
title('Sample Trace: TFEB-EYFP');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 4 3]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}