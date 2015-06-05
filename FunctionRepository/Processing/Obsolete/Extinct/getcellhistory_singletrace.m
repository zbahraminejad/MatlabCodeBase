%%%%%%%% Review images and save single cell traces %%%%%%%%%%%%%%%%%%%%
% Author: Mingyu Chung
% Last Revision: 9/4/2012
%
% Use:
% 1. Run calcnuclei_saveimages to save cropped & processed images
% 2. Run tracenuclei
% 
% This code requires: 
% - subaxis.m
% - datatip_image_singletrace_angie.m
% - datatip_image_singletrace_ctd.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;
global path moviename allx ally allxorg allyorg plotfignum

%%%%%% directory paths %%%%%%%%%%%%%%%%%%
path = 'h:\Documents\Timescape\20120807_12drugs\';
datadir = [path,'Data\'];
savedir = [path,'GoodCells\'];
%%%%%% User settings %%%%%%%%%%%%%%%%%%%%
selectmode = 1;         %1=enable imaging and selection. 0=viewing traces only
sample = [];            %can specify a set of cells.  set empty vector to consider all cells
start = [];             %enter starting cellid (empty defaults to 1)
ymin = 0.3; ymax = 2; yrange = ymax-ymin;
historycutoff = 0.7;
%%%%%%%%%%%% Experiment Settings %%%%%%%%
moviename = '3_9_1';
numframes = 208;
framesperhr = 5;
drugtime = 90;          %set to zero if no drug
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([datadir,moviename,'_alldata'],'bestsp','best_rc')
%% store signal values in angieratio and cdt1 matrices
totalcells = size(bestsp{end},1);
totalframes = size(bestsp,3);
angieratio = -10000*ones(totalcells,totalframes);
cdt1 = -10000*ones(totalcells,totalframes);
for f = 1:size(angieratio,2)
    tempcell = find(bestsp{f}(:,1)~=0);                                     %ignore any cells that have negative coordinates (I'm not sure how negative coordinates would ever appear anyways...)
    angieratio(tempcell,f) = bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);  %cytosol/nuclear ratio for DHB (Angie's) sensor
    cdt1(tempcell,f) = bestsp{f}(tempcell,6);                               %median nuclear Cdt1 value
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
for i = 1:totalcells
    %%% inherit all ancestor cell history and note each mitosis
    temptrace = i;
    m = 1;
    while best_rc(temptrace,5)~=best_rc(temptrace,2)
        mother = best_rc(temptrace,2);
        daughterfirst = best_rc(temptrace,1);
        motherfirst = best_rc(mother,1);
        best_rc_temp(i,1) = motherfirst;         %reset cell's first frame to mother, but use temp until finished
        goodmitoses(i,m) = daughterfirst;
        m = m+1;
        angieratio(i,daughterfirst) = mean([angieratio(i,daughterfirst) angieratio(mother,daughterfirst)]);
        angieratio(i,motherfirst:daughterfirst-1) = angieratio(mother,motherfirst:daughterfirst-1);
        for f=motherfirst:daughterfirst-1
            allx(i,f) = bestsp{f}(mother,1);
            ally(i,f) = bestsp{f}(mother,2);
        end
        temptrace = mother;
    end
    %%% find each direct descendent mitosis
    daughtermitoses = best_rc(find(best_rc(:,2)==i & best_rc(:,1)<=best_rc(i,3)),1);  %ignore cells that split (falsely) after end of cell history
    if m==1         %best_rc col 2 = col 5 for first mother, so remove
        daughtermitoses(find(daughtermitoses==best_rc(i,1))) = [];
    end
    goodmitoses(i,m:m+length(daughtermitoses)-1) = daughtermitoses;
end
best_rc_org = best_rc;
best_rc = best_rc_temp;

%% plotting
set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);
counter = 0;
j = 0;
tc_sel = [1:totalcells];
if sample
    tc_sel = sample;
end
savetracelist = zeros(totalcells,1);      %will store all selected cell trace id#s
if start
    tc_sel = [start:totalcells];
    load([savedir,moviename,'_goodcells_temp.mat'],'tracelist','apoptosislist');
    savetracelist(1:length(tracelist)) = tracelist;
end
for counter=1:length(tc_sel)
    cellid = tc_sel(counter);
    if best_rc(cellid,3)-best_rc(cellid,1)+1 >= totalframes*historycutoff  % fraction of movie that at least one member of the group must exist for 
        if selectmode
            clf;
            set(gcf,'Position',[round(0.02*screenx) round(0.05*screeny) round(0.95*screenx) round(0.4*screeny)]);
            plotfignum = gcf;           %reference this in datatip_image
        else
            if sample
                subaxis(1,length(tc_sel),counter, 'SH', 0.04, 'SV', 0.1,'Margin',0.05);
            else
                j = j+1;
                figure(ceil(j/24));             %only 24 plots per figure
                subplot(4,6,mod(j-1,24)+1);
            end
        end
        set(gcf,'color','w');
        set(gca,'YTick',0.5:0.5:2,'XTick',0:framesperhr*4:numframes,'Position',[0.03 0.1 0.95 0.8]);
        axis([1,numframes,ymin,ymax]);
            
        angieratioeachcell = angieratio(cellid,:);     %get angieratio for each frame for this cell
        cdt1eachcell = cdt1(cellid,:);                 %same with cdt1
        if best_rc(cellid,1)~=best_rc_org(cellid,1)
            x1 = best_rc(cellid,1):best_rc_org(cellid,1)-1;    %ancestor history
            y1 = angieratioeachcell(x1);
            x2 = best_rc_org(cellid,1):best_rc(cellid,3);      %this cell's history
            y2 = angieratioeachcell(x2);
            line(x1,y1,'LineStyle','--','DisplayName',num2str(cellid));
            line(x2,y2,'DisplayName',num2str(cellid));
            x = [x1,x2];
            y = [y1,y2];
        else
            x = best_rc(cellid,1):best_rc(cellid,3);                %x is start frame through end frame for cell cc
            y = angieratioeachcell(x);                              %y contains intensity thru time of fluorescent protein you wish to plot
            traceangie = line(x,y,'DisplayName',num2str(cellid));   %this plots the YFP intensity through time 
        end
        %%%%% mark all mitoses
        mitoses = goodmitoses(cellid,:);
        mitoses = mitoses - best_rc(cellid,1)+1;
        mitoses = mitoses(find(mitoses>0));
        for cc = 1:length(mitoses)
            hold on;
            plot(x(mitoses(cc)),y(mitoses(cc)), 'go', 'markerfacecolor', 'g')   %plot a dot at frame of automatically found mitosis
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
                %set(dcm_obj,'Updatefcn',@datatip_traceid);          %don't require special fcn anymore
                response = input('','s');
                figure(gcf);
                if strcmp(response,'y')      %save traceid: 'y','enter'
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
        save([savedir,moviename,'_goodcells_temp.mat'],'tracelist');
    end   
end

