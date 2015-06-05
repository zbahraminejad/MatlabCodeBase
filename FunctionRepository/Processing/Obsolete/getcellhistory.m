%function getcellhistory(row,col,site)
%%%%%%%% Review images and save single cell traces %%%%%%%%%%%%%%%%%%%%
% Author: Mingyu Chung
% Last Revision: 9/4/2012
%
% Use:
% Before running, execute the following code:
% 1. calcnuclei_saveimages: save cropped & processed images, and store
% features for each cell in each image
% 2. tracenuclei: track cells
% 
% This code requires:
% - subaxis.m
% - datatip_image_singletrace_angie.m
% - datatip_image_singletrace_ctd.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
%clear; close all;
global path moviename allx ally allxorg allyorg plotfignum

%%%%%% directory paths %%%%%%%%%%%%%%%%%%
path = 'h:\Documents\Timescape\20120807_12drugs\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120202_6p5-30p5hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120209_24-48hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120115_48-72hr\';
datadir = [path,'Data\'];
savedir = [path,'GoodCells\'];
nocall = 1;
if nocall   %running without parallelwell.m (i.e. not calling fxn)
    row = 4;
    col = 11;
    site = 0;
end
%%%%%%%%%%%% Experiment Settings %%%%%%%%
framesperhr = 5;
drugtime = 90;          %set to zero if no drug
%%%%%% User settings %%%%%%%%%%%%%%%%%%%%
selectmode = 1;         %1=enable imaging and selection. 0=viewing traces only
sample = [60];          %[4,11,0] oscillations in hDHB & Cdt1
%sample = [109];         %[4,11,1]
%selectin = [2 46 48 96 217 383 418 436 447 449 465 468 470 502 505 510 518 525 539 550 552 561 683]; %[4 11 1] errant mitosis
selectin = [];
start = [];             %enter starting cellid (empty defaults to 1) in case re-starting manual selection
directgating = 0;       %1=skip plots and save all good traces. 0=standard mode
registry = 1;
selectout = [];
plotwindows = 0;
plotsignal2 = 1;
plotsignal2peak = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin = 0; ymax = 2.5; yrange = ymax-ymin;
historycutoff = 0.5;
screening = 1;
dotsize = 8;
signalwidth = 2;        %linewidth of normal trace
windowlinewidth = 2;    %linewidth over specified windows of time
timewindowmatrix = [1,20;21,40;41,60;61,80];
palette = 'cymk';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

moviename = [num2str(row),'_',num2str(col),'_',num2str(site)];
if registry
    load([datadir,moviename,'_alldata'],'bestsp','best_rc','corruptlist','leaveoutlist')
else
    load([datadir,moviename,'_alldata'],'bestsp','best_rc')
end
%% store signal values in angieratio and cdt1 matrices
totalcells = size(bestsp{end},1);
totalframes = size(bestsp,3);
cellidlist = [1:totalcells];
if registry
    leavein = cellidlist(~ismember(cellidlist,leaveoutlist));                   %enables tracking cell ids through tracenuclei pgm
else
    leavein = cellidlist;
end
signal1 = -10000*ones(totalcells,totalframes);
signal2 = -10000*ones(totalcells,totalframes);
nucsig = -10000*ones(totalcells,totalframes);
angieratio = -10000*ones(totalcells,totalframes);
for f = 1:size(signal1,2)
    tempcell = find(bestsp{f}(:,1)~=0);                                     %ignore any cells that have negative coordinates (I'm not sure how negative coordinates would ever appear anyways...)
    signal1(tempcell,f) = bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);     %cytosol/nuclear ratio for DHB (Angie's) sensor
    %signal1(tempcell,f) = bestsp{f}(tempcell,3);
    signal2(tempcell,f) = bestsp{f}(tempcell,6)/300;                            %median nuclear Cdt1 value
    %signal2(tempcell,f) = bestsp{f}(tempcell,4)/200;
    nucsig(tempcell,f) = bestsp{f}(tempcell,3);
    %nucsiz(tempcell,f) = bestsp{f}(tempcell,4);
    angieratio(tempcell,f) = bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);
end
best_rc_temp = best_rc;
%%%% complete all traces and record mitoses
allx = zeros(totalcells,totalframes);
ally = allx;
allxorg = allx;
allyorg = allx;
for f=1:totalframes
    allxorg(1:size(bestsp{f},1),f) = bestsp{f}(:,1);
    allyorg(1:size(bestsp{f},1),f) = bestsp{f}(:,2);
end
allx = allxorg;
ally = allyorg;
goodmitoses = zeros(totalcells,10);
for i = leavein
    %%% inherit all ancestor cell history and note each mitosis
    temptrace = i;
    m = 1;
    while best_rc(temptrace,5)~=best_rc(temptrace,2)
        mother = best_rc(temptrace,2);
        daughterfirst = best_rc(temptrace,1);
        motherfirst = best_rc(mother,1);
        best_rc_temp(i,1) = motherfirst;         %reset cell's first frame to mother, but use temp until finished
        mothermitoses = best_rc(leavein(find(best_rc(leavein,2)==mother & best_rc(leavein,1)<daughterfirst)),1); %record all prior mitoses of mother
        mitnum = length(mothermitoses);
        goodmitoses(i,m:m+mitnum-1) = mothermitoses;
        m = m+mitnum;
        signal1(i,motherfirst:daughterfirst-1) = signal1(mother,motherfirst:daughterfirst-1);
        signal2(i,motherfirst:daughterfirst-1) = signal2(mother,motherfirst:daughterfirst-1);
        nucsig(i,motherfirst:daughterfirst-1) = nucsig(mother,motherfirst:daughterfirst-1);
        angieratio(i,motherfirst:daughterfirst-1) = angieratio(mother,motherfirst:daughterfirst-1);
        for f=motherfirst:daughterfirst-1
            allx(i,f) = bestsp{f}(mother,1);
            ally(i,f) = bestsp{f}(mother,2);
        end
        temptrace = mother;
    end
    %%% find each direct descendent mitosis
    daughtermitoses = best_rc(leavein(find(best_rc(leavein,2)==i)),1);
    daughtermitoses = unique(daughtermitoses);
    daughtermitoses(find(daughtermitoses==best_rc(i,1))) = [];
    goodmitoses(i,m:m+length(daughtermitoses)-1) = daughtermitoses;
    if m>1
        goodmitoses(i,find(goodmitoses(i,:)==best_rc(mother,1))) = 0;
    end
end
best_rc_org = best_rc;
best_rc = best_rc_temp;

%% screen out cells with negative angie values
for cc=1:length(leavein)
    i = leavein(cc);
    if min(angieratio(i,best_rc(i,1):best_rc(i,3)))<0
        leavein(cc) = 0;
    end
end
leavein = leavein(leavein>0);

%% screen out cells with missed or errant mitoses
mitosiserror = zeros(totalcells,1);
scanwidth = 7;
%leaveinuncorrupted = leavein(~ismember(leavein,corruptlist));
for i = leavein
    %%%% infer where mitoses should be by nuc signal, nuc size, and angie
    tracklength = [best_rc(i,1):best_rc(i,3)];
    nucsigdrop = [];
    for tt = 1:scanwidth
        tracklengthmod = tracklength(mod(tracklength,scanwidth)+1==tt);
        nucsigdrop = [nucsigdrop,tracklengthmod(diff(nucsig(i,tracklengthmod))<-0.3)]; %-0.3
    end
    nucsigdrop = unique(nucsigdrop);
    nucsigdropzone = zeros(1,totalframes);
    if isempty(nucsigdrop)
        nucsigdrop = [];
    end
    for k = nucsigdrop
        window = scanwidth-1;
        priorsig = 0;
        if k>5
            priorsig = nucsig(i,k-5);
        end
        windowsignal = nucsig(i,k:k+window);
        if windowsignal(1)==max(windowsignal) & windowsignal(1)>priorsig & windowsignal(end)<mean(windowsignal([1 end]))
            nucsigdropzone(k:k+window) = 1;
        end
    end
    
    mitosissignatures = nucsigdropzone;
    mitoses = goodmitoses(i,:);
    mitoses = mitoses(mitoses>0);
    for j = mitoses
        if mitosissignatures(j)==0
            mitosiserror(i) = 1;
        else
            mitosissignatures(j)=0;
            checkright = j; checkleft = j;
            if j<length(mitosissignatures)
                checkright = j+1;
            end
            if j>1
                checkleft = j-1;
            end
            while mitosissignatures(checkright)==1 | mitosissignatures(checkleft)==1
                if mitosissignatures(checkright)==1
                    mitosissignatures(checkright)=0;
                    if checkright<length(mitosissignatures)
                        checkright = checkright+1;
                    end
                end
                if mitosissignatures(checkleft)==1
                    mitosissignatures(checkleft)=0;
                    if checkleft>1
                        checkleft = checkleft-1;
                    end
                end
            end
        end
    end
    if max(mitosissignatures)>0
        mitosiserror(i) = 1;
    end
end
mitosiserrorlist = cellidlist(find(mitosiserror));

%% plotting
set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);
counter = 0;
j = 0;
tc_sel = leavein;
if sample
    tc_sel = sample;
end
savetracelist = zeros(1,totalcells);      %will store all selected cell trace id#s
if start
    tc_sel = tc_sel(find(tc_sel==start):end);
    load([savedir,moviename,'_goodcells.mat'],'tracelist');
    savetracelist(1:length(tracelist)) = tracelist;
end
for counter=1:length(tc_sel)
    cellid = tc_sel(counter);
    %%%%%% Gate on cell history %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    enoughhistory = best_rc(cellid,3)-best_rc(cellid,1)+1 >= totalframes*historycutoff;
    %enoughhistory = best_rc(cellid,1)<50 & best_rc(cellid,3)==totalframes;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if registry
        corrupted = ismember(cellid,corruptlist);
    end
    errantmitosis = ismember(cellid,mitosiserrorlist);
    selectedout = ismember(cellid,selectout);
    selectedin = ismember(cellid,selectin);
    selection = 1;
    if screening
        %selection = enoughhistory & ~errantmitosis & ~corrupted & ~selectedout;
        selection = enoughhistory & errantmitosis & ~selectedout;
        %selection = enoughhistory;
    end
    if selection
        j=j+1;
        if directgating
            nextindex = find(~savetracelist,1);
            savetracelist(nextindex) = cellid;
            continue                    %skip all plotting and just record gated cell
        end
        if selectmode
            clf;
            %set(gcf,'Position',[round(0.02*screenx) round(0.05*screeny) round(0.95*screenx) round(0.4*screeny)]);
            set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
            plotfignum = gcf;           %reference this in datatip_image
            set(gca,'Position',[0.03 0.1 0.95 0.8]);
        else
            if sample
                subaxis(1,length(tc_sel),counter, 'SH', 0.04, 'SV', 0.1,'Margin',0.05);
            else
                figure(ceil(j/24));             %only 24 plots per figure
                subplot(4,6,mod(j-1,24)+1);
            end
        end
        set(gcf,'color','w');
        set(gca,'YTick',0.5:0.5:2,'XTick',0:framesperhr*4:totalframes);
        axis([1,totalframes,ymin,ymax]);
        signal1eachcell = signal1(cellid,:);                    %get angieratio for each frame for this cell
        signal2eachcell = signal2(cellid,:);                    %same with cdt1
        if best_rc(cellid,1)~=best_rc_org(cellid,1)
            x1 = best_rc(cellid,1):best_rc_org(cellid,1)-1;     %ancestor history
            y1 = signal1eachcell(x1);
            x2 = best_rc_org(cellid,1):best_rc(cellid,3);       %this cell's history
            y2 = signal1eachcell(x2);
            line(x1,y1,'LineStyle','--','DisplayName',num2str(cellid));
            line(x2,y2,'DisplayName',num2str(cellid));
            x = [x1,x2];
            ysig1 = [y1,y2];
        else
            x = best_rc(cellid,1):best_rc(cellid,3);            %x is start frame through end frame for cell cc
            ysig1 = signal1eachcell(x);                         %y contains intensity thru time of fluorescent protein you wish to plot
            ysig1 = smooth(ysig1);
            line(x,ysig1,'DisplayName',num2str(cellid),'linewidth',signalwidth);        %this plots the YFP intensity through time 
        end
        
        
        %%%%% mark all mitoses %%%%%%%%%%%%%%
        mitoses = goodmitoses(cellid,:);
        mitoses = mitoses(mitoses>=best_rc(cellid,1));
        relativemitoses = mitoses - best_rc(cellid,1)+1;
        if length(relativemitoses)>0
            for cc = 1:length(relativemitoses)
                hold on;
                if relativemitoses(cc)<=length(x)
                    plot(x(relativemitoses(cc)),ysig1(relativemitoses(cc)), 'ro', 'markerfacecolor', 'r','markersize',dotsize)   %plot a dot at frame of automatically found mitosis
                end
            end
        end
        if plotsignal2
            ysig2 = signal2eachcell(x);
            ysig2 = smooth(ysig2);
            line(x,ysig2,'Linestyle','--','color','g','linewidth',signalwidth);
        end
        firstmitosisindex = find(mitoses>=drugtime,1);
        if firstmitosisindex
            firstmitosis = mitoses(firstmitosisindex);
            framemax = x(end);
            if firstmitosisindex < length(mitoses)
                framemax = mitoses(firstmitosisindex+1);
            end
            firstmitosisrel = firstmitosis - best_rc(cellid,1)+1;
            framemax = framemax - best_rc(cellid,1)+1;
            
            if plotwindows
                for tw = 1:size(timewindowmatrix,1)
                    timewindow = firstmitosisrel+timewindowmatrix(tw,1):firstmitosisrel+timewindowmatrix(tw,2);
                    twend = max(timewindow(timewindow<=framemax));
                    timewindow = timewindowmatrix(tw,1):timewindowmatrix(tw,2);
                    twend = max(timewindow(timewindow<=framemax-firstmitosisrel));
                    if twend
                        line(x(firstmitosisrel+timewindow(1):firstmitosisrel+twend),ysig1(firstmitosisrel+timewindow(1):firstmitosisrel+twend),'color',palette(tw),'linewidth',windowlinewidth);
                    end
                end
            end
            if plotsignal2peak
                [maxval,maxindex] = max(ysig2(firstmitosisrel:framemax));
                hold on;
                if maxindex < framemax-firstmitosisrel+1-3
                    plot(x(maxindex+firstmitosisrel-1),maxval,'ko','markerfacecolor','k','markersize',dotsize,'linewidth',signalwidth);
                end
            end
        end
     
        title(['Cell ',num2str(cellid)]);
        if drugtime
            line([drugtime drugtime],[ymin ymax],'Color','r');  %add vertical line when drug is added
        end
        %%%% selectmode: view and choose cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 'enter': see next cell group
        % click on a trace: 
        %       'y','enter': saves cell id
        %       'a','enter': saves cell id & adds to apoptosis list
        % 'i','enter': enter image viewing mode
        %   * must 'delete' previous datatips prior to entering image mode.
        %   - click on a trace to see image. arrow keys update next point.
        % 'd','enter': return to data viewing mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if selectmode
            dcm_obj = datacursormode(gcf);
            datacursormode on;
            response = 'on';
            while response                                          %get next graph: 'enter'
                response = input('','s');
                figure(gcf);
                if strcmp(response,'y')                             %save traceid: 'y','enter'
                    nextindex = find(~savetracelist,1);
                    savetracelist(nextindex) = cellid;
                    fprintf(['saving cell ',num2str(cellid),'\n']);
                end
                if strcmp(response,'a')                             %view angie sensor: 'a','enter'
                    delete(findall(gcf,'Type','hggroup'));          %must remove datapoints before changing program
                    set(dcm_obj,'Updatefcn',@datatip_image_singletrace_angie);
                    while ~strcmp(response,'d')                     %exit image mode:'d','enter'
                        response = input('','s');
                    end
                    delete(findall(gcf,'Type','hggroup'));
                    close(gcf+1);                                   %close the image window
                end
                if strcmp(response,'c')                             %view ctd sensor: 'c','enter'
                    delete(findall(gcf,'Type','hggroup'));          %must remove datapoints before changing program
                    set(dcm_obj,'Updatefcn',@datatip_image_singletrace_ctd);
                    while ~strcmp(response,'d')                     %exit image mode:'d','enter'
                        response = input('','s');
                    end
                    delete(findall(gcf,'Type','hggroup'));
                    close(gcf+1);                                   %close the image window
                end
                figure(gcf);
            end
        end
    end
    if selectmode
        tracelist = savetracelist(find(savetracelist));
        save([savedir,moviename,'_goodcells.mat'],'tracelist');
    end
end
fprintf([num2str(j),'\n']);
if directgating
    tracelist = savetracelist(find(savetracelist));
    save([savedir,moviename,'_goodcells.mat'],'tracelist');
end
if directgating | selectmode
    save([savedir,moviename,'_singletracedata.mat'],'signal1','signal2','goodmitoses');
end
%set(gcf,'color','w','PaperPosition',[0 0 16 6]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');
cd ..\Processing; %return to this directory