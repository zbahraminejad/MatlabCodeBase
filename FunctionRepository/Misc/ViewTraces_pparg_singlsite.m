%  function ViewTraces
row=2;col=2;site=1;

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
% projectpath='D:\Documents\Projects\';
imagepath='G:\Michael\';
%imagepath='E:\';
experimentpath='20150401-CC-Diff\';
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
%     rawdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=8;
framesperhr=5.5;
drugspike=0/framesperhr;
drugspike2 = 0/framesperhr;
startframe = 1;
endframe = 653;
frames=startframe:endframe;
channelnames={'Cy5_' 'YFP_' 'mCherry_'};
edgemask=1; %0:no masks saved 1: masks saved
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 0;    %0:one trace per plotn 1:all traces in one plot
selectmode = 1;  %View trace images
selectin=[960];   %Choose trace IDs (necessary for selectmode)
plotsignal2 = 1;
plotsignal3 = 0;
IFoption=0; %0:no IFdata 1:IFdata
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1 =200000;
ymin2=0; ymax2=3; %ymax2=1000;
%%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=0; %0:no gating 1:daug1hters 2:no daughters
quiescentanalysis=0;
%%%%%% Load Data %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datafile =[datadir,'tracedata_',shot,'_nolink',a','genealogy','jitters');'.mat'];
datafile =[datadir,'tracedata_',shot,'_nolink','.mat'];
% load([datadir,'tracedata_',shot,'.mat'],'tracedat
load(datafile,'tracedata','genealogy','jitters');
% [tracedata,tracestats,motherstats,IFdata,IDs]=gathertracedata_3(datadir,datafile,shot,motheroption,daughteroption,IFoption);
[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother]=gathertracedata_mz_1(datadir,datafile,shot,motheroption,daughteroption,IFoption);
% hist(tracedata(:,end,7),50)
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
% reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
%%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[tempcellcount,~]=size(tracestats);
alteredlength  = zeros(tempcellcount,1);
alteredstart = zeros(tempcellcount,1);
gatelength = alteredlength;
for i = 1:tempcellcount
    alteredlength(i) = sum(~isnan(tracedata(i,:,1)));
    alteredstart(i) = find(~isnan(tracedata(i,:,1)),1,'first');
    gatelength(i) = sum(~isnan(tracedata(i,99:end,1)));
end
tracestats(:,1) = alteredstart;
tracestats(:,3) = alteredlength;
minlengthtrace=300;
badlengths=gatelength<minlengthtrace;
%%% gate PPARg data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelPPARg=10;nucareachannel = 3;
noisethresh = 10^10;
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=10^10;  %threshold above which max of each trace must be %50
minpos = 0;      %0:anywhere 1:firstframe 2:lastframe
minthresh = 0; %threshold below which min of each trace must be %50
% [traces1,badtraces1]=gate_pparg_2_singlesite(tracedata,tracestats,noisethresh,channelPPARg,nucareachannel,minthresh,minpos,maxthresh,maxpos);
[traces1,badtraces1]=gate_pparg_1_allsites(tracedata,tracestats,noisethresh,channelPPARg,nucareachannel,minthresh,minpos,maxthresh,maxpos);
% [traces2,badtraces2]=gate_pparg_8_mother(tracedata,tracestats,noisethresh,channelPPARg,minthresh,minpos,maxthresh,maxpos);
%%% gate Degron data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelGem=11; 
maxpos=0;     %0:anywhere 1:firstframe 2:lastframe
maxthresh=100;  %threshold above which max of each trace must be %50
minpos=0;      %0:anywhere 1:mothertrace 2:daughtertrace
minthresh=25; %threshold below which min of each trace must be %50
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
% orphan = isnan(tracestats(:,4));
badstart = tracestats(:,1)> 100;
overlap = zeros(tempcellcount,1);
repeatedcelltrace = unique(lastcellmother,'sorted');
repeatedcelltrace(repeatedcelltrace==-1)=[];
overlap(repeatedcelltrace) = 1; overlap(lastcellmother==-1)=0;
%%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces = badlengths | badtraces1 | badstart | overlap;% | badtraces2;% | orphan;  
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
% traces3=traces3(~badtraces,:);
tracedata=tracedata(~badtraces,:,:);
tracestats=tracestats(~badtraces,:);
markedmitosis = markedmitosis(~badtraces);
IDs=IDs(~badtraces,:);

%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
% traces1=normalizetraces_3(traces1,tracestats,1);
traces2=normalizetraces_3(traces2,tracestats,1);
% traces3=normalizetraces_3(traces3,tracestats,0);

numgated=size(tracestats,1);
% if numgated>192
%     numgated=96;
% end
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
steps=5; ystep1=round(((ymax1-ymin1)/steps)*10)/10; ystep2=round(((ymax2-ymin2)/steps)*10)/10;
trace2color=[0 0 1];
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
%     drugtime2=drugspike2;
%     drugplot2=0;
else
    drugtime=0;
    drugplot=drugspike;
%     drugtime2=0;
%     drugplot2=drugspike2;
end
xtime=frames/framesperhr;
% xtime=xtime-drugtime;
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
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(traces1(i,:),6); %default 3
    if ~plotsignal2
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
    elseif plotsignal2
        ysig2=smoothignorenans(traces2(i,:),5);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1,'YColor','k','FontName','Arial','FontSize',8);
        set(hline1,'DisplayName',num2str(i),'color',[0.9 .75 0],'linewidth',tracewidth);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike>0
        line([drugplot drugplot],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
        %line([drugtime drugtime],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
    end
    
    if drugspike2>0
%         line([drugplot2 drugplot2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
        line([drugtime2 drugtime2],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','-');
    end
    %%%%% mark mitosis (only the current cell's birth) %%%%%%%%%%%%%%%%%%%%
    if ~isnan(tracestats(i,4))
        plot(xtime(markedmitosis{i}),ysig1(markedmitosis{i}),'ro','markerfacecolor', 'g','markersize',dotsize);
    end
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2,'FontName','Arial','FontSize',8);
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if plotsignal3
        ysig3=smoothignorenans(traces3(i,:),3);
        line(xtime,ysig3,'color','g','DisplayName',num2str(i),'linewidth',tracewidth);
    end
    title(strcat('IDs=',num2str(IDs(i)),', row=', num2str(row),', col=',num2str(col),', site=',num2str(site)));
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