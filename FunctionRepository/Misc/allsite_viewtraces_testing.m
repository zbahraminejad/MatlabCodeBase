row=2;col=5;
numsites = 3;
alltraces1 = [];
alltracestats = [];
alltraces2 = [];
allIDs = [];
for site = 1:numsites

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
% projectpath='D:\Documents\Projects\';
imagepath = 'D:\Michael\';
%imagepath='E:\';
experimentpath='12072014-Michael-48hr-CC\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
wellname = nameandsite(shot);
%shot=wellnum2str(row,col,site);
datadir=([imagepath,experimentpath,'Data\']);
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=16;
framesperhr=4;
drugspike=0/framesperhr;
frames=1:187;
channelnames={'CFP_' 'YFP_' 'mCherry_'};
edgemask=1; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble=0;    %0:one trace per plotn 1:all traces in one plot
selectmode=0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2=0;
plotsignal3=0;
IFoption=0; %0:no IFdata 1:IFdata
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1=1.5;
ymin2=0; ymax2=1.5; %ymax2=1000;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=1; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
[tracedata,tracestats,motherstats,IFdata,IDs]=gathertracedata_3(datadir,shot,motheroption,daughteroption,IFoption);
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
% reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace=40;
badlengths=tracestats(:,3)<minlengthtrace;
%%% gate degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelGem=7;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=40;  %threshold above which max of each trace must be %50
minpos=0;      %0:anywhere 1:mothertrace 2:daughtertrace
minthresh=20; %threshold below which min of each trace must be %50
[traces1,badtraces1]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
% [traces2,badtraces2]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
%%% gate degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelCdt1=6;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=100;  %threshold above which max of each trace must be %50
minpos=0;      %0:anywhere 1:mothertrace 2:daughtertrace
minthresh=20; %threshold below which min of each trace must be %50
% [traces1,badtraces1]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelCdt1,maxthresh,minthresh,maxpos,minpos);
[traces2,badtraces2]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelCdt1,maxthresh,minthresh,maxpos,minpos);
%%% gate miscellaneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%traces3=tracedata(:,:,3); %Area

%%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badstarttime = tracestats(:,1)>150;
badtraces= badlengths | badtraces1; %| badtraces2; % | badstarttime;
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
%traces3=traces3(~badtraces,:);
tracedata=tracedata(~badtraces,:,:);
tracestats=tracestats(~badtraces,:);
IDs=IDs(~badtraces,:);
% cclength = tracestats(:,3)./4
%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
traces1=normalizetraces_3(traces1,tracestats,0);
traces2=normalizetraces_3(traces2,tracestats,0);
%traces3=normalizetraces_3(traces3,tracestats,0);

alltraces1 = [alltraces1; traces1];
alltraces2 = [alltraces2;traces2];
alltracestats = [alltracestats;tracestats];
allIDs = [allIDs; IDs];
end

%%%%%%%%for cdt1 analysis%%%%%%%%
% [h,w]=size(alltraces1);
% numtraces = h;
% alignedtraces = zeros(h,w);
% firstnonzeroind = zeros(numtraces,1);
% toobig = zeros(numtraces,1);
% for i = 1:numtraces
%     
%     alignpoint = alltracestats(i,1);
%     alignend = 187;
%     alignedtraces(i,1:length(alignpoint:alignend)) = alltraces1(i,alignpoint:alignend);
%     alignedtraces(i,length(alignpoint:alignend)+1:187)=nan;
%     lastnonzeroind = find(~isnan(alignedtraces(i,:)),1,'last');
%     if alignedtraces(i,lastnonzeroind) > 0.1
%         toobig(i) = 1;
%     else
%         toobig(i) = 0;
%     end
% end
% alltracestats(logical(toobig),:)=[];
% alignedtraces(logical(toobig),:) = [];
% alltracestats(alignedtraces(:,1)>0.5,:)=[];
% alignedtraces(alignedtraces(:,1)>0.5,:)=[];
% alltraces1 = alignedtraces;
% [numtraces,~] = size(alltraces1);
% traceduration = zeros(numtraces,1);
% for trace = 1:numtraces
%    sig = alltraces1(trace,:);
%    smoothsig = smoothignorenans(sig,15);
%    [~,maxindex] = max(smoothsig);
%    temptrace = smoothsig(maxindex:end);
%    [~,minloc] = min(temptrace);
%    finalloc = minloc+maxindex;
%    traceduration(trace) = finalloc/4;
% end
% traceduration

%%%%%%%%%%for geminin analysis%%%%%%%%%%%%%%%%%%
[numtraces,~] = size(alltraces1);
badpks = zeros(numtraces,1);
duration = zeros(numtraces,1);

for trace = 1:numtraces
sig = alltraces1(trace,:);
smoothsig = smoothignorenans(sig,8); % smooth trace
smoothsig(smoothsig<0) = 0;
[allpks,allpklocs] = findpeaks(smoothsig,'minpeakdistance',20); % find all peaks and their corresponding location
realpks = allpks(allpks > 0.2); % threshold peaks based on intensity
realpklocs = allpklocs(allpks>0.2);
numpks = length(realpks);
%break signal up into numpeaks+1 parts and find the minimum
switch numpks
    case 0 
        badpks(trace) = 1;
        duration(trace) = -1;
    case 1
        [min1,minloc1] = min(smoothsig(1:realpklocs)); 
%         [min2,minloc2] = min(smoothsig(realpklocs+1:end)); add1 = length(1:realpklocs); %need to correct for distance scale)
        if min1 < 0.1 %&& min2 < 0.1
         badpks(trace) = 0;
         duration(trace) = (realpklocs - minloc1)/framesperhr;
        else
         badpks(trace) = 1;
         duration(trace) = -1;
        end
    case 2 
%         [min1,minloc1] = min(smoothsig(1:realpklocs(1)));
%         [min2,minloc2] = min(smoothsig(realpklocs(1)+1:realpklocs(2))); add1 = length(1:realpklocs(1));
%         [min3,minloc3] = min(smoothsig(realpklocs(2)+1:end)); add2 = add1+length(realpklocs(1)+1:realpklocs(2));
%         allmins = [min1 min2 min3];
%         allminlocs = [minloc1 minloc2+add1 minloc3+add2];
%         allminlocs(allmins>0.1) = [];
%         diffmin = diff(allminlocs);
%         if isempty(diffmin)
            badpks(trace) = 1;
            duration(trace) = -1;
%         else            
%             badpks(trace) = 0;
%             duration(trace) = diffmin(1)/framesperhr;
%         end
    case 3
%         [min1,minloc1] = min(smoothsig(1:realpklocs(1)));
%         [min2,minloc2] = min(smoothsig(realpklocs(1)+1:realpklocs(2))); add1 = length(1:realpklocs(1));
%         [min3,minloc3] = min(smoothsig(realpklocs(2)+1:realpklocs(3))); add2 = add1+length((realpklocs(1)+1):realpklocs(2));
%         [min4,minloc4] = min(smoothsig(realpklocs(3)+1:end)); add3 = add2 + length((realpklocs(2)+1):realpklocs(3));
%         allmins = [min1 min2 min3 min4];
%         allminlocs = [minloc1 minloc2+add1 minloc3+add2 minloc4+add3];
%         allminlocs(allmins>0.1) = [];
%         diffmins = diff(allminlocs);
%         if isempty(diffmins)
            badpks(trace) = 1;
            duration(trace) = -1;
%         else
%             badpks(trace) = 0;
%             duration(trace) = diffmins(1)/framesperhr;
%         end
end
end
% duration(logical(badpks)) = [];
% 
% alltracestats(logical(badpks),:)=[];
% alltraces1(logical(badpks),:)=[];
%%%%%%%%%%%%%% PLOTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
numgated=size(alltracestats,1);
% if numgated>192
%     numgated=192;
% end
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
steps=5; ystep1=round((ymax1-ymin1)/steps); ystep2=round((ymax2-ymin2)/steps);
trace2color=[1 0 0];
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
else
    drugtime=0;
    drugplot=drugspike;
end
xtime=frames/framesperhr;
xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=find(ismember(IDs,selectin));
end
if ensemble
    figure; hold on;
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
%         figure(site)
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(alltraces1(i,:),10); %default 3
    
%     ysig1 = [ysig1(1),ysig1];
%     ysig1 = diff(ysig1);
    if ~plotsignal2
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
    elseif plotsignal2
        ysig2=smoothignorenans(alltraces2(i,:),10);
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
    if ~isnan(alltracestats(i,4))
        plot(xtime(alltracestats(i,1)),ysig1(alltracestats(i,1)),'ro','markerfacecolor', 'g','markersize',dotsize);
    end
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2);
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if plotsignal3
        ysig3 = smoothignorenans(traces3(i,:),3);
        line(xtime,ysig3,'color','g','DisplayName',num2str(i),'linewidth',tracewidth);
    end
    title(num2str(allIDs(i)));
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