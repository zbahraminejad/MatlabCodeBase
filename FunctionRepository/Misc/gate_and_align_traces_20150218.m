% clear all 
% close all
rowmat = [2,3,4;5,6,7];
rows= rowmat(1,:); col=4;
% rows = 7;
sites = 1:2;
numrows = length(rows);
numsites = length(sites);
alltraces1 = [];
alltracestats = [];
alltraces2 = [];
allIDs = [];    
allmotherstats = [];
allstaindata = [];
allLength = [];
for row = 1:numrows
    for site = 1:numsites

    global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
    % projectpath='D:\Documents\Projects\';
    imagepath = 'D:\Michael\';
    %imagepath='E:\';
    % experimentpath='11202014-Michael-CellCycle-48hour-glass\';
    experimentpath='12052014-CellCycle-24hr-48hr-bins\';
    % experimentpath='102420140-CellCycle-48hr\';
    shot=[num2str(rows(row)),'_', num2str(col), '_', num2str(sites(site))];
    wellname = nameandsite(shot);
    %shot=wellnum2str(row,col,site);
    datadir=([imagepath,experimentpath,'Data\']);
    separatedirectories=0;
    if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
    else
    %     rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
    rawdir = [imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    end
    immunoframe=0;

    %%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucr=12;
    framesperhr=6;
    drugspike=0/framesperhr;
    frames=2:139; %319
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
    motheroption=0; %0:no gating 1:mothers (traces that end in mitosis) 2:no mothers (traces that don't end in mitosis)
    daughteroption=0; %0:no gating 1:daughters (traces that start with mitosis) 2:no daughters (traces that don't start with mitosis)
    quiescentanalysis=0;
    removequiescent = 0; % 0 for no gating, 1 to look for cells that stay quiescent the whole time;
    %%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    datafile = [datadir,'tracedata_',shot,'.mat'];
    load(datafile,'tracedata','genealogy','jitters');
    [tracedata,tracestats,motherstats,IFdata,IDs]=gathertracedata_3(datadir,datafile,shot,motheroption,daughteroption,IFoption);
    tracedata(:,1,:) = []; %this line and the mislabed experiment
    tracestats(:,1:2)=tracestats(:,1:2)-1; % Subtract one for the 24-48 hour bin (-1) ;
    %%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
    % reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
    %%% gate length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     lengthgating = tracedata(:,:,1);
%     [numcells,~] = size(lengthgating);
%     totallength = zeros(numcells,1);
%     for cell = 1:numcells
%        totallength(cell) = sum(~isnan(lengthgating(cell,:))); 
%     end
    minlengthtrace=60;
    badlengths=tracestats(:,3)<minlengthtrace;%| isnan(tracedata(:,1,1));
%     badlengths = totallength<minlengthtrace;
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
    noisethresh=0.2;  %threshold for max positive jump in cyto/nuc ratio
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
    [numcells,~ ]=size(traces1);
%     badstart = traces1(:,1) > 0.35 & tracestats(:,1)==1;
%     badstart = isnan(traces1(:,6));
    %%% combine all gates and remove %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     disp(length(tracestats)-sum(badlengths))
    badtraces= badlengths | badtraces1;%|badtraces2;% | badIF | orphan;%| badtraces2;% | badstart;%  badIF | orphan; %badtraces1 | badtraces2;% | badtraces2;
    traces1=traces1(~badtraces,:);
    traces2=traces2(~badtraces,:);
    % traces3=traces3(~badtraces,:);
    tracedata=tracedata(~badtraces,:,:);
    tracestats=tracestats(~badtraces,:);
    IDs=IDs(~badtraces,:);
%     staindata = staindata(~badtraces);

    %%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %0:normalize each trace by max of its own trace
    %1:normalize by median of max of all traces if min>max*0.5
    %2:normalize each trace by value at first frame after mitosis
    % traces1=normalizetraces_3(traces1,tracestats,1);
    traces2=normalizetraces_3(traces2,tracestats,0);
    % traces3=normalizetraces_3(traces3,tracestats,0);

    alltraces1 = [alltraces1; traces1];
    alltraces2 = [alltraces2;traces2];
    alltracestats = [alltracestats;tracestats];
    allIDs = [allIDs; IDs];
    allmotherstats = [allmotherstats; motherstats];
    [cellcount,~] = size(traces1);
    allLength = [allLength cellcount];
%     allstaindata = [allstaindata; staindata];
    end
 
end


[numcells,~] = size(alltraces1);
g1start = zeros(numcells,1);
spstart = zeros(numcells,1);
wrongtrace = zeros(numcells,1);
for cell = 1:numcells
   start = alltracestats(cell,1); % start of the trace or start of cell life after anaphase (depends on cell)
   smoothed = smoothignorenans(alltraces1(cell,:),6);
   diffsmoothed = diff(smoothed);
   maxval = max(smoothed(start+12:end));
   badhigh = smoothed(start+5)>0.4; %& start<6 | smoothed(start)>1.5 & start >=6;
%    badlow = maxval<1.3 ; 
   badslope = median(diffsmoothed(1:6)) > 0.01;
   if  badhigh || maxval <1.2
       g1start(cell) = start;
       wrongtrace(cell) = 1;
       continue
   end
   tracefrag = smoothed(start:end);
   [~,minInd] = min(tracefrag(1:round(length(tracefrag)/2)));
   offset = start - 1;
   maxInd = find(tracefrag(minInd:end)>=0.8,1,'first');
   maxInd = minInd + (maxInd-1);
   if any([isempty(maxInd) isempty(minInd) isempty(tracefrag)])
       g1start(cell) = start;
       wrongtrace(cell) = 1;
       continue
   end
   tracefrag = tracefrag(1:maxInd); % trace we will work with temporarily
   fragdiff = [tracefrag(1) tracefrag];
   tempdiff = diff(fragdiff,1,2);
   tempdiff = smooth(tempdiff,6);
%    if tempdiff(1)>0.01
%        g1start(cell) = 1+offset+start;
%    else
%    hitindices = find(tempdiff>0.007,20,'first');
   hitindices = find(tempdiff>0.015);
   hitindices(hitindices<minInd) = [];
   hitindices(hitindices>length(smoothed)-12) = [];
   if numel(hitindices) > 1
      tempstore = zeros(numel(hitindices),1);
      for hit = 1:numel(hitindices)
          absindices = hitindices+offset;
          temptrace = diff(smoothed(absindices(hit)+1:absindices(hit)+10));
          highscore = temptrace; %- smoothed(absindices(hit));
          highscore(highscore>0.01) = 1;
          highscore(highscore<0.01&highscore>0) = 0;
          highscore(highscore<0) = 0;
          tempstore(hit) = sum(highscore);
      end
      tempind = find(tempstore>=4,1,'first');
      finalindex = hitindices(tempind);
   end
%    limitind = find(smoothtemp>min(smoothtemp) + (max(smoothtemp)-min(smoothtemp))/2,1,'first');
%    hitindices(hitindices>limitind)=[];
%    difflist = abs(tempdiff(hitindices) - .01);
%    [~,ind] = min(difflist);
   g1start(cell) = finalindex+offset;  
%    end
%    gemtrace = alltraces2(cell,g1start(cell):end);
%    if gemtrace(1) > 0.01
%        wrongtrace(cell) = 1;
%        continue
%    end
%    smoothgem = smoothignorenans(gemtrace,6);
%    cut = length(smoothgem) - 10;
%    smoothgem = [smoothgem smoothgem(end)];
%    diffgem = diff(smoothgem,1,2);
%    gemhit = find(diffgem>0.01);
%    gemhit(gemhit>cut) =[];
%    if numel(gemhit)>1
%        tempstore = zeros(numel(gemhit),1);
%        for hit = 1:numel(gemhit)
%             tempstore(hit) = mean(smoothgem(gemhit(hit):gemhit(hit)+5));        
%        end
%        ind = find(tempstore>0.075,1,'first');
%        gemhit = gemhit(ind);
%    end
%    spstart(cell) = g1start(cell) + gemhit-1;
end
wrongtrace = logical(wrongtrace);


% figure,hist(cclength,20)

% [time,alltraces1] = aligntraces_5(alltraces1,alltracestats(:,1),alltracestats,allmotherstats,daughteroption);
% [time,alltraces1] = aligntraces_5(alltraces1,g1start,alltracestats,allmotherstats,daughteroption);
% badcdk2start = alltraces1(:,time==0);
% badcdk2start = badcdk2start>0.6;

% [time,alltraces2] = aligntraces_5(alltraces2,g1start,alltracestats,allmotherstats,daughteroption);
[numcells,~] = size(alltraces1);
% colorcondition = zeros(numcells,1);
% colorcondition(1:sum(allLength(1:6))) = 1;
% colorcondition = logical(colorcondition);

%%%%%%Search for quiescent cells%%%%%%%%%%%%%%
% hightemp = zeros(numcells,1);
% quiescent = zeros(numcells,1);
% for i = 1:numcells
%     cycling = sum(alltraces1(i,alltracestats(i,1):end)>1);
%     maxval = max(alltraces1(i,alltracestats(i,1):end));
%     tooshort = alltracestats(i,3)<120;
%     if cycling >10 || maxval > 1 %|| tooshort
%         quiescent(i) = 0;
%     else
%         quiescent(i) = 1;
%     end
%     hightemp(i) = any(alltraces1(i,alltracestats(i,1):alltracestats(i,1)+5)>0.4);
% end
%%%%%%Last minute Gating%%%%%%%%%%%
secondgate = wrongtrace;% | badcdk2start;logical(hightemp);
alltracestats(secondgate,:) = [];
alltraces1(secondgate,:) = [];
alltraces2(secondgate,:) = [];
g1start(secondgate) = [];
spstart(secondgate) = [];
% quiescent(secondgate) = [];


%%%% cell cycle timing %%%%%%%%%%
cclength = (alltracestats(:,2)-g1start)./6;
g1length = (spstart-g1start)./6;
intermitotic = (g1start-alltracestats(:,1))./6;
[mean(cclength) median(cclength) std(cclength)]
% [mean(g1length) median(g1length) std(g1length)]
% [mean(intermitotic) median(intermitotic) std(intermitotic)]
% [sum(quiescent)+length(g1length) sum(quiescent) length(g1length)]

testsort = sort(g1start);
[testsort,testind] = sort(g1start);
sortedtraces = alltraces1(testind,:);
% figure,imagesc(sortedtraces)
% colorbar
%% 
% imagesc(alltraces1, [0 3])
numgated=size(alltraces1,1);
if numgated>192
    numgated=48;
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
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
%         figure(ceil(counter/6));
%         subaxis(2,3,mod(counter-1,6)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.04); %5x4  ML, MR, MT, MB control outter borders and SH controls spacing between subplots
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(alltraces1(i,:),12); %default 3
%     ysig1 = [ysig1(1) ysig1];
%     ysig1 = diff(ysig1);
%     ysig1 = smoothignorenans(ysig1,12);
    if ~plotsignal2
%         if quiescent(counter) == 1
%             linecolor = [1 0 0];
% %         elseif quiescent(counter) <0
% %             linecolor = [0 1 0];
%         else
%             linecolor = [0 0 1];
%         end
        linecolor = [0 0 1];
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
    plot(xtime(g1start(i)),ysig1(g1start(i)),'ko','markerfacecolor', 'k','markersize',dotsize);
%     plot(xtime(spstart(i)),ysig2(spstart(i)),'bo','markerfacecolor', 'b','markersize',dotsize);
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