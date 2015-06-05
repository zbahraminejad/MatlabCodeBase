%--------------------------------------------------------------------------
% Function ViewTraces
%----------------------------------NOTES-----------------------------------
% 'tracedata' represents a three dimensional data structure with the
% following format: [cell_identity, tracing_frames, imaging_channels]
%--------------------------------------------------------------------------
% 'tracestats' represents a two dimensional data structure with the
% following format: row: cell_identity, columns: start_frame, end_frame,
% trace_duration, genealogy
%--------------------------------------------------------------------------
% CDK2 Data is represented by 'traces1' data structure, obtained using
% gate_Cdk2_8_mother(...) function.
%--------------------------------------------------------------------------
% Degron Data is represented by 'traces2' data structure, obtained using
% gate_Geminin_8_mother(...) function.
%--------------------------------------------------------------------------
% Furthermore, 'traces3' data structure is equivalent to tracedata(:,:,4)
% which represents total h2b intensity.
%--------------------------------------------------------------------------

close all;
clearvars

% Specify the row, column, and site corresponding to
% the experiment that trace viewer is going to display
row = 6; col = 9; site = 1;

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr

%% File Paths Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagePath = 'G:\Zahra\';
experimentPath = '150507_LiveCell\';
shot = [num2str(row), '_', num2str(col), '_', num2str(site)];
wellName = nameandsite(shot);
% To be used when data is supplied through a local drive
% dataDir = ([imagePath,experimentPath,'Data\']);
% To be used when data is supplied through a USB removable drive
dataDir = '/Volumes/KINGSTON/Data/';
separateDirectories = 0;

if separateDirectories == 1
    rawdir = [imagePath,experimentPath,'Raw\',shot,'\']; % Path to Raw data files
	maskdir = [imagePath,experimentPath,'Mask\',shot,'\']; % Path to mask data files
else
	rawdir = [imagePath,experimentPath,'Raw\',wellName,shot,'_']; % Path to Raw data files
    maskdir = [imagePath,experimentPath,'Mask\',shot,'_']; % Path to mask data files
end

immunoframe = 0;

%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr = 8;
framesperhr = 6;
drugSpikeOne = 0/framesperhr;
drugSpikeTwo = 0/framesperhr;
startFrame = 1;
endFrame = 561;
frames = startFrame:endFrame;
channelnames = {'Cy5_' 'YFP_' 'mCherry_'};
edgeMask = 1; % Definition: '0' means no masks saved, '1' means masks saved

%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ensemble = 0;         % 0: one trace per plot 1: all traces in one plot
selectMode = 0;       % View trace images
selectin = [];        % Choose trace IDs (necessary for selectMode)
plotSignalTwo = 0;
plotSignalThree = 0;
IFoption = 0;         % 0: no IFdata 1: IFdata

%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1 = 0; ymax1 = 80000;
ymin2 = 0; ymax2 = 3;

%% Cell-Cycle Analysis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motherOption = 0;      % 0: no gating 1: mothers (traces that end in mitosis) 2: no mothers (traces that don't end in mitosis)
daughterOption = 0;    % 0: no gating 1: daughters (traces that start with mitosis) 2: no daughters (traces that don't start with mitosis)
quiescentAnalysis = 0;
onlyQuiescent = 0;     % 0: no gating, 1: to look for cells that stay quiescent the whole time (in a state or period of inactivity or dormancy)
removequiescent = 0;

%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFile = [dataDir,'tracedata_',shot,'_nolink','.mat'];
load(dataFile,'tracedata','genealogy','jitters');

% Extract gated trace data, trace statistics, and mother statistics
% [tracedata,tracestats,motherStats,IFdata,IDs] = gathertracedata_3(dataDir,dataFile,shot,motherOption,daughterOption,IFoption);
[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother] = gathertracedata_mz_1(dataDir,dataFile,shot,motherOption,daughterOption,IFoption);
% hist(tracedata(:,end,7),50)
% tracedata(:,1,:) = []; % this line and the mislabed experiment
% tracestats(:,1:2) = tracestats(:,1:2)-1; % Subtract one for the 24-48 hour bin (-1);

%% For Purpose of Visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tracestats = getstats(tracedata,genealogy); %start,end,length,genealogy
% reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);

%% Gate Length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minLengthTrace = 400;

[tempCellCount,~]=size(tracestats);
alteredLength  = zeros(tempCellCount,1);
alteredStart = zeros(tempCellCount,1);
for i = 1:tempCellCount
    alteredLength(i) = sum(~isnan(tracedata(i,:,1)));
    alteredStart(i) = find(~isnan(tracedata(i,:,1)),1,'first');
end
tracestats(:,1) = alteredStart;
tracestats(:,3) = alteredLength;

badLengths = tracestats(:,3) < minLengthTrace;

%% Smooth, post-Process, and Gate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% PPARg Data
%--------------------------------------------------------------------------
% Depending on the experiment, this needs to be updated
% channelPPARg = 8; % Michael's experiments
channelPPARg = 6; % Zahra's experiments
nucareaChannel = 3;
noiseThresh = 10^10;
maxpos = 0;         % 0: anywhere 1: first frame 2: last frame
maxThresh = 10^10;  % threshold above which maximum of each trace must be % 50
minpos = 0;         % 0: anywhere 1: mother trace 2: daughter trace
minThresh = 0;      % threshold below which minimum of each trace must be % 50
% [traces1,badtraces1] = gate_pparg_2_singlesite(tracedata,tracestats,noiseThresh,channelPPARg,nucareaChannel,minthresh,minpos,maxthresh,maxpos);
[traces1,badtraces1] = gate_pparg_1_allsites(tracedata,tracestats,noiseThresh,channelPPARg,nucareaChannel,minThresh,minpos,maxThresh,maxpos);

%--------------------------------------------------------------------------
% Degron Data
%--------------------------------------------------------------------------
% channelGem = 10;
% maxpos = 0;      % 0: anywhere 1: first frame 2: last frame
% maxThresh = 40;  % threshold above which maximum of each trace must be % 50
% minpos = 0;      % 0: anywhere 1: mother trace 2: daughter trace
% minThresh = 20;  % threshold below which minimum of each trace must be % 50
% [traces1,badtraces1] = gate_Geminin_8_mother(tracedata,tracestats,motherStats,channelGem,maxThresh,minThresh,maxpos,minpos);
% tic
% [traces2,badtraces2] = gate_Geminin_8_mother(tracedata,tracestats,motherStats,channelGem,maxThresh,minThresh,maxpos,minpos);
% toc

%--------------------------------------------------------------------------
% CDK2 Data
%--------------------------------------------------------------------------
% nucChannel = 6; cytoChannel = 8;
% nucThreshOption = 0;
% 0: normalize each trace by percentile of entire trace
% 1: normalize by first frame of trace (regardless of whether mitosis or not)
% nucThresh = 20;     % threshold for nuc intensity according to nucThreshOption

% motherThresh = 0;   % threshold for max DHB ratio of mother. default gating = 1.0. 0: ignore
% noiseThresh = 0.2;  % threshold for max positive jump in cyto/nuc ratio
% tic
% [traces1,badtraces1] = gate_Cdk2_8_mother(tracedata,nucChannel,cytoChannel,tracestats,nucThreshOption,nucThresh,motherThresh,noiseThresh,quiescentAnalysis,onlyQuiescent);
% toc
% [traces2,badtraces2] = gate_Cdk2_8_mother(tracedata,nucChannel,cytoChannel,tracestats,nucThreshOption,nucThresh,motherThresh,noiseThresh,quiescentAnalysis,onlyQuiescent);

%% Gate by IF Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% staindata = IFdata(:,8);
% thresholdsig = 300;
% badIF = (staindata < thresholdsig);
% orphan = isnan(tracestats(:,4)) & tracestats(:,3) < 100;
% figure(6),hist(staindata,20)

%% Gate Miscellaneous %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% traces3 = tracedata(:,:,4); % Area is 3, Mass(total h2b intensity) is 4
% traces2 = tracedata(:,:,7);
% traces3 = tracedata(:,:,8);

%% Combine All Gates and Remove Garbage Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces = badLengths | badtraces1; %| badtraces2; %| badIF | orphan; % badtraces1 | badtraces2; %| badtraces2;
traces1 = traces1(~badtraces,:);
% traces2 = traces2(~badtraces,:);
% traces3 = traces3(~badtraces,:);
tracedata = tracedata(~badtraces,:,:);
tracestats = tracestats(~badtraces,:);
IDs = IDs(~badtraces,:);

%% Normalize Traces by Maximum in Lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0: normalize each trace by maximum of its own trace
% 1: normalize by median of maximum of all traces if minimum > maximum*0.5
% 2: normalize each trace by value at first frame after mitosis
% traces1=normalizetraces_3(traces1,tracestats,1);
% traces2=normalizetraces_3(traces2,tracestats,1);
% traces3=normalizetraces_3(traces3,tracestats,0);

numGated = size(tracestats,1);
if numGated > 192
    numGated = 192;
end

%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims = get(0,'ScreenSize'); screenx = screendims(3); screeny = screendims(4);

%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
windowSize = 20;        % Moving average window size
movingAverageType = ['Simple','Exponential','Triangular','Weighted','Modified'];
dotsize = 8;            % default = 8 presentation = 12
mainTraceWidth = 2;     % default = 2 presentation = 4
estimateTraceWidth = 1; % default = 0.5 presentation = 1
steps = 5; ystep1 = round(((ymax1-ymin1)/steps)*10)/10; ystep2 = round(((ymax2-ymin2)/steps)*10)/10;
trace2color = [1 0 0];

if selectMode == 0
    drugTimeOne = drugSpikeOne;
    % drugTimeTwo = drugSpikeTwo;
    drugPlotOne = 0;
    % drugPlotTwo = 0;
else
    drugTimeOne = 0;
    % drugTimeTwo = 0;
    drugPlotOne = drugSpikeOne;
    % drugPlotTwo = drugSpikeTwo;
end

xtime = frames/framesperhr;
% xtime = xtime - drugTimeOne;
selection = 1:numGated;

if ~isempty(selectin)
    selection = find(ismember(IDs,selectin));
end

%% Reduce Memory Footprint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars tracedata motherstats genealogy frames IFdata jitters lastcellmother
clearvars alteredLength alteredStart badLengths badtraces badtraces1

%% Structure Holding Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ySignalStruct = [];

%% Cell-Cycle Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ensemble % Draw and hold all traces in just one plot
    figure; hold on;
end

% Iterate over all cell traces
for counter = 1:length(selection)
    i = selection(counter);
    if selectMode
        % Clear current figure window
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        % Get current figure handle
        plotfignum = gcf;
        % set(gca,'Position',[0.03 0.1 0.95 0.8]);
    elseif ~ensemble
        figure(ceil(counter/24));
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
        % figure(ceil(counter / 4));
        % subaxis(2,2,mod(counter-1,4)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.05); %5x4
    end

    set(gcf,'color','w');
    ysig1 = smoothignorenans(traces1(i,:),3); % default: 3

    if ~plotSignalTwo
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',mainTraceWidth);
        line(xtime,tsmovavg(smoothignorenans(ysig1,13)','s',windowSize,1),...
            'DisplayName',num2str(i),'linewidth',estimateTraceWidth, 'Color', 'red');
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        % Store the output signal for future use
        ySignalStruct = setfield(ySignalStruct,strcat('ysig1_',num2str(IDs(i)),...
            '_', num2str(row),'_',num2str(col),'_',num2str(site)),ysig1);
        % Plot 2-D line plots with y-axes on both left and right side
        % Calculate and plot the simple moving average
        % [haxes,hline1,hline2] = plotyy(xtime,ysig1,xtime,tsmovavg(smoothignorenans(ysig1,13)','s',windowSize,1));
        % axes(haxes(1));
        % axis([xtime(1) xtime(end) ymin1 ymax1]);
        % Calculate and plot the exponential weighted moving average moving average
        % [haxes,hline1,hline2] = plotyy(xtime,ysig1,xtime,tsmovavg(ysig1','e',windowSize,1));
		% Calculate and plot the triangular moving average moving average
		% [haxes,hline1,hline2] = plotyy(xtime,ysig1,xtime,tsmovavg(ysig1','t',windowSize,1));
		% Calculate and plot the weighted moving average moving average
		% semiGaussian = [0.026 0.045 0.071 0.1 0.12 0.138];
		% semiGaussian = [semiGaussian fliplr(semiGaussian)];
        % [haxes,hline1,hline2] = plotyy(xtime,ysig1,xtime,tsmovavg(ysig1','w',semiGaussian,1));
		% Calculate and plot the modified moving average moving average
		% [haxes,hline1,hline2] = plotyy(xtime,ysig1,xtime,tsmovavg(ysig1','m',windowSize,1));
    elseif plotSignalTwo % Draw the graph together in one plot
        ysig2 = smoothignorenans(traces2(i,:),3);
		% Plot 2-D line plots with y-axes on both left and right side
        [haxes,hline1,hline2] = plotyy(xtime,ysig1,xtime,ysig2);
        ySignalStruct = setfield(ySignalStruct,strcat('ysig2_',num2str(IDs(i)),...
            '_', num2str(row),'_',num2str(col),'_',num2str(site)),ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        % set(hline1,'DisplayName',num2str(i),'color','b','linewidth',mainTraceWidth);
        set(hline1,'DisplayName',num2str(i),'color',[0.9 .75 0],'linewidth',mainTraceWidth);
    end

    hold on;

    %%%%% Mark Drug Spike and Title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugSpikeOne > 0
        line([drugPlotOne drugPlotOne],[ymin1 ymax1],'Color','k','linewidth',mainTraceWidth,'linestyle','--');
        % line([drugTimeOne drugTimeOne],[ymin1 ymax1],'Color','k','linewidth',mainTraceWidth,'linestyle','--');
    end

    if drugSpikeTwo > 0
        % line([drugPlotTwo drugPlotTwo],[ymin1 ymax1],'Color','k','linewidth',mainTraceWidth,'linestyle','-');
        line([drugTimeTwo drugTimeTwo],[ymin1 ymax1],'Color','k','linewidth',mainTraceWidth,'linestyle','-');
    end
    %%%%% Mark Mitosis (only the current cell's birth) %%%%%%%%%%%%%%%%%%%%
    if ~isnan(tracestats(i,4))
		% Start frame of a given cell's trace determines the instance mitosis has occurred
        plot(xtime(markedmitosis{i}),ysig1(markedmitosis{i}),'ro','markerfacecolor', 'g','markersize',dotsize);
        % plot(xtime(tracestats(i,1)),ysig1(tracestats(i,1)),'ro','markerfacecolor', 'g','markersize',dotsize);
    end

    %%% Adjust Signal2 Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotSignalTwo
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2);
        % set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',mainTraceWidth);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2,'FontName','Arial','FontSize',8);
    end

	%%% Adjust Signal3 Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotSignalThree
        % ysig3 = smoothignorenans(traces3(i,:),3);
        % line(xtime,ysig3,'color','g','DisplayName',num2str(i),'linewidth',mainTraceWidth);
        trace3color = [0 0.8 0];
        ysig3 = smoothignorenans(traces3(i,:),3);
        ySignalStruct = setfield(ySignalStruct,strcat('ysig3_',num2str(IDs(i)),...
            '_', num2str(row),'_',num2str(col),'_',num2str(site)),ysig3);
        line(xtime,ysig3,'color',trace3color,'DisplayName',num2str(i),'linewidth',mainTraceWidth);
    end
    % title(num2str(IDs(i)));
    title(strcat('ID: ',num2str(IDs(i)),', Row: ',num2str(row),', Col: ',num2str(col),', site: ',num2str(site)));

    %%% Operate Image-Viewing Mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectMode
        selectmodedisplay(edgeMask);
    end
end

%% Moving Average %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% % Iterate over all cell traces for moving averages
% for counter = 1:length(selection)
%     i = selection(counter);
%     set(gcf,'color','w');
% 
%     ysig1 = smoothignorenans(traces1(i,:),3); % default: 3
%     if plotSignalTwo
%         ysig2 = smoothignorenans(traces2(i,:),3);
%     end
% 
% 	% Calculate and plot variety of moving averages for the supplied data
% 	if ~ensemble
% 		[sigOneHeight,sigOneWidth] = size(ysig1);
% 		sigOneMovingAverage = zeros(sigOneHeight,sigOneWidth,length(movingAverageType));
% 		% Calculate the simple moving average
% 		sigOneMovingAverage(1,:,movingAverageType(1)) = tsmovavg(ysig1','s',windowSize,1);
% 		% Calculate the exponential weighted moving average moving average
% 		sigOneMovingAverage(1,:,movingAverageType(2)) = tsmovavg(ysig1','e',windowSize,1);
% 		% Calculate the triangular moving average moving average
% 		sigOneMovingAverage(1,:,movingAverageType(3)) = tsmovavg(ysig1','t',windowSize,1);
% 		% Calculate the weighted moving average moving average
% 		semiGaussian = [0.026 0.045 0.071 0.1 0.12 0.138];
% 		semiGaussian = [semiGaussian fliplr(semiGaussian)];
% 		sigOneMovingAverage(1,:,movingAverageType(4)) = tsmovavg(ysig1','w',semiGaussian,1);
% 		% Calculate the modified moving average moving average
% 		sigOneMovingAverage(1,:,movingAverageType(5)) = tsmovavg(ysig1','m',windowSize,1);
% 		% Plot the results for the five moving average calculations
% 		plot(xtime,ysig1,...
% 			xtime,sigOneMovingAverage(1,:,movingAverageType(1)),...
% 			xtime,sigOneMovingAverage(1,:,movingAverageType(2)),...
% 			xtime,sigOneMovingAverage(1,:,movingAverageType(3)),...
% 			xtime,sigOneMovingAverage(1,:,movingAverageType(4)),...
% 			xtime,sigOneMovingAverage(1,:,movingAverageType(5)))
% 		legend('CDK Data','Simple','Exponential','Triangular','Weighted',...
% 			'Modified','Location','NorthWest')
% 		title('CDK Data under Different Moving Averages Regimes')
% 
%
% 		if plotSignalTwo
% 			[sigTwoHeight,sigTwoWidth] = size(ysig2);
% 			sigTwoMovingAverage = zeros(sigTwoHeight,sigTwoWidth,length(movingAverageType));
% 			% Calculate the simple moving average
% 			sigTwoMovingAverage(1,:,movingAverageType(1)) = tsmovavg(ysig1','s',windowSize,1);
% 			% Calculate the exponential weighted moving average moving
% 			% average
% 			sigTwoMovingAverage(1,:,movingAverageType(2)) = tsmovavg(ysig1','e',windowSize,1);
% 			% Calculate the triangular moving average moving average
% 			sigTwoMovingAverage(1,:,movingAverageType(3)) = tsmovavg(ysig1','t',windowSize,1);
% 			% Calculate the weighted moving average moving average
% 			semiGaussian = [0.026 0.045 0.071 0.1 0.12 0.138];
% 			semiGaussian = [semiGaussian fliplr(semiGaussian)];
% 			sigTwoMovingAverage(1,:,movingAverageType(4)) = tsmovavg(ysig1','w',semiGaussian,1);
% 			% Calculate the modified moving average moving average
% 			sigTwoMovingAverage(1,:,movingAverageType(5)) = tsmovavg(ysig1','m',windowSize,1);
% 			% Plot the results for the five moving average calculations
% 			plot(xtime,ysig1,...
% 				xtime,sigTwoMovingAverage(1,:,movingAverageType(1)),...
% 				xtime,sigTwoMovingAverage(1,:,movingAverageType(2)),...
% 				xtime,sigTwoMovingAverage(1,:,movingAverageType(3)),...
% 				xtime,sigTwoMovingAverage(1,:,movingAverageType(4)),...
% 				xtime,sigTwoMovingAverage(1,:,movingAverageType(5)))
% 				legend('Degron Data','Simple','Exponential','Triangular','Weighted',...
% 					'Modified','Location','NorthWest')
% 				title('Degron Data under Different Moving Averages Regimes')
% 		end
%
% 	end
% 
% 	hold on;
% 
%     %%%%% Mark Drug Spike and Title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if drugSpikeOne > 0
%         % line([drugPlotOne drugPlotOne],[ymin1
%         % ymax1],'Color','k','linewidth',mainTraceWidth,'linestyle','--');
%         line([drugTimeOne drugTimeOne],[ymin1 ymax1],'Color','k','linewidth',mainTraceWidth,'linestyle','--');
%     end
% end

fprintf([num2str(numGated),'\n']);

%
%  title('Sample Trace: TFEB-EYFP');
%  xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
%  set(gcf,'color','w','PaperPosition',[0 0 4 3]); %3x4 or mini
%  saveas(gcf,'h:\Downloads\Fig.jpg');
%  close(gcf);

