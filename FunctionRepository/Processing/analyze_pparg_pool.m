rows=3:5;cols=6;
numsites = 4;
numreps = length(rows);
reps = rows;
alltraces1 = [];
alltracestats = [];
alltraces2 = [];
allIDs = [];
% row=5;col=7;site=1;
global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
for rep = 1:numreps
for site = 1:numsites
% projectpath='D:\Documents\Projects\';
imagepath='E:\Michael\';
%imagepath='E:\';
experimentpath='20150126-gem-pparg-diff\';
shot=[num2str(reps(rep)),'_', num2str(cols), '_', num2str(site)]; % reps(rep) can either be row or column depending on experiment setup
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
%     rawdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=8;
framesperhr=6;
drugspike=0/framesperhr;
drugspike2 = 0/framesperhr;
startframe = 1;
endframe = 584;
frames=startframe:endframe;
channelnames={'mcherry_' 'cfp_' 'yfp_'};
edgemask=1; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 0;    %0:one trace per plotn 1:all traces in one plot
selectmode = 0;  %View trace images
selectin=[];   %Choose trace IDs (necessary for selectmode)
plotsignal2=1;
plotsignal3=0;
IFoption=0; %0:no IFdata 1:IFdata
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1 =800;
ymin2=0; ymax2=500; %ymax2=1000;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafile =[datadir,'tracedata_',shot,'_nolink','.mat'];
% load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
load(datafile,'tracedata','genealogy','jitters');
[tracedata,tracestats,motherstats,IFdata,IDs]=gathertracedata_3(datadir,datafile,shot,motheroption,daughteroption,IFoption);
% hist(tracedata(:,end,7),50)
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
% reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minlengthtrace=100;
badlengths=tracestats(:,3)<minlengthtrace;
%%% gate degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelPPARg=7;
noisethresh = 100;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=10000;  %threshold above which max of each trace must be %50
minpos = 0;      %0:anywhere 1:firstframe 2:lastframe
minthresh = 0; %threshold below which min of each trace must be %50
[traces1,badtraces1]=gate_pparg_8_mother(tracedata,tracestats,noisethresh,channelPPARg,minthresh,minpos,maxthresh,maxpos);
% [traces2,badtraces2]=gate_pparg_8_mother(tracedata,tracestats,noisethresh,channelPPARg,minthresh,minpos,maxthresh,maxpos);
%%% gate degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelGem=6;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=100;  %threshold above which max of each trace must be %50
minpos=0;      %0:anywhere 1:mothertrace 2:daughtertrace
minthresh=50; %threshold below which min of each trace must be %50
% [traces1,badtraces1]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
[traces2,badtraces2]=gate_Geminin_8_mother(tracedata,tracestats,motherstats,channelGem,maxthresh,minthresh,maxpos,minpos);
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nucchannel=7; cytochannel=9;
% nucthreshoption=0;
% %0:normalize each trace by percentile of entire trace
% %1:normalize by first frame of trace (regardless of whether mitosis or not)
% nucthresh=100;    %threshold for nuc intensity according to nucthreshoption
% 
% motherthresh=0;   %threshold for max DHB ratio of mother. default gating=1.0.  0:ignore
% noisethresh=0.2;  %threshold for max positive jump in cyto/nuc ratio
% [traces2,badtraces2]=gate_Cdk2_8_mother(tracedata,nucchannel,cytochannel,tracestats,nucthreshoption,nucthresh,motherthresh,noisethresh,quiescentanalysis);
%%% gate miscellaneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% traces3=tracedata(:,:,4); % Area is 3, Mass(total h2b intensity) is 4
% traces2 = tracedata(:,:,7);
% traces3 = tracedata(:,:,8);
%%% gate out cells that didn't have a last frame around day 4%%%%%%%%%%%%%%
endcells = tracedata(:,endframe-18:endframe,1);
total = sum(~isnan(endcells),2);
badend = total==0;
orphan = isnan(tracestats(:,4)) & tracestats(:,3)<250;

%%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces=badlengths | badend | orphan; %|badtraces1 | badtraces2;
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
% traces3=traces3(~badtraces,:);
tracedata=tracedata(~badtraces,:,:);
tracestats=tracestats(~badtraces,:);
IDs=IDs(~badtraces,:);

%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
% traces1=normalizetraces_3(traces1,tracestats,1);
% traces2=normalizetraces_3(traces2,tracestats,1);
% traces3=normalizetraces_3(traces3,tracestats,0);
alltraces1 = [alltraces1; traces1];
alltraces2 = [alltraces2;traces2];
alltracestats = [alltracestats;tracestats];
allIDs = [allIDs; IDs];
end
end

numgated=size(alltracestats,1);
% numgated=size(tracestats,1);
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