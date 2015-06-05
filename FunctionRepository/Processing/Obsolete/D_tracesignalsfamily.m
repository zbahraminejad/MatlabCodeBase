clear; close all;
selectmode = 0;
global path moviename bestsp best_rc traceid plotfignum
path = 'h:\Documents\Timescape\20120807_12drugs\';
datadir = [path,'Data\'];
savedir = [path,'GoodCells\'];
%%%%%%%%%%%% Experiment Settings %%%%%%%%
moviename = '3_9_1';
numframes = 208;
framesperhr = 5;
drugtime = 90;          %set to zero if no drug
sample = [];          %set empty vector to show all groups
%%%%%%%%%%%% Other Settings %%%%%%%%%
ymin = 0.3; ymax = 2; yrange = ymax-ymin;
historycutoff = 0.7;

load([datadir,moviename,'_alldata'],'bestsp','best_rc')
%% store signal values in angieratio and cdt1 matrices
totalcells = size(bestsp{end},1);
angieratio = -10000*ones(totalcells,size(bestsp,3));
cdt1 = -10000*ones(totalcells,size(bestsp,3));
for f = 1:size(angieratio,2)
    tempcell = find(bestsp{f}(:,1)~=0);                                     %ignore any cells that have negative coordinates (I'm not sure how negative coordinates would ever appear anyways...)
    angieratio(tempcell,f) = bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);  %cytosol/nuclear ratio for DHB (Angie's) sensor
    cdt1(tempcell,f) = bestsp{f}(tempcell,6);                               %median nuclear Cdt1 value
end
%remove negative values and assign them to be 0
%angieratio(angieratio<0 & angieratio~=-10000)=0;
%cdt1(cdt1<0 & cdt1~=-10000)=0;

%% grouping
GP = zeros(totalcells,1);          %no difference between size(bestsp,1) and size(best_rc,1)...
counter = 0;
for cc = 1:size(GP,1)
    if best_rc(cc,2) ~= best_rc(cc,5)   %split occurred (from mother to daughter)
        GP(cc) = GP(best_rc(cc,2));     %give the daughter the same GP value as mother (which occurred earlier, so already has a GP# assigned)
    else
        counter = counter+1;
        GP(cc) = counter;
    end
end

%% plotting
%iptsetpref('ImshowBorder','tight');     %removes borders around figures
set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);
%set(0,'Units','normalized');            %normalizes screen size for future positioning
counter = 0;
colors = 'ygbkr';
tc_sel = [1:max(GP)];
if sample
    tc_sel = sample;
end
savetracelist = zeros(totalcells,1);      %will store all selected cell trace id#s

j=0;
for counter=1:length(tc_sel)  
    gp_rc = best_rc(GP==tc_sel(counter),:);                 %store best_rc of all row(cell)#s in GP that have value 'tc'
    gp_rc(:,6) = mod(1:size(gp_rc,1),5)+1;                  %to loop through the 7 colors and give each cell a color in the 6th col of gp_rc0
    if (max(gp_rc(:,3))-min(gp_rc(:,1))+1) >= (size(bestsp,3)*historycutoff)  % fraction of movie that at least one member of the group must exist for 
        if selectmode
            clf;
            %set(gcf,'Position',[round(0.7*screenx) round(0.6*screeny) round(0.28*screenx) round(0.3*screeny)]);       %place plot in upper left
            set(gcf,'Position',[round(0.02*screenx) round(0.05*screeny) round(0.95*screenx) round(0.4*screeny)]);
            plotfignum = gcf;           %reference this in datatip_image
            set(gca,'Position',[0.03 0.1 0.95 0.8]);
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
        set(gca,'YTick',0.5:0.5:2,'XTick',0:framesperhr*4:numframes);
        axis([1,numframes,ymin,ymax]);
        setx = zeros(size(gp_rc,1),1);
        sety = setx;
        for cc = 1:size(gp_rc,1)                           %for each cell      
            cellid = gp_rc(cc,5);
            tempcolor = colors(gp_rc(cc,6));
            angieratioeachcell = angieratio(cellid,:);     %get angieratio for each frame for this cell
            cdt1eachcell = cdt1(cellid,:);                 %same with cdt1
            x = gp_rc(cc,1):gp_rc(cc,3);                   %x is start frame through end frame for cell cc
            y = angieratioeachcell(x);                     %y contains intensity thru time of fluorescent protein you wish to plot
            traceangie = line(x,y,'color',tempcolor,'DisplayName',num2str(cellid));  %this plots the YFP intensity through time 
            if cc > 5
                set(traceangie,'LineStyle','-.');
            end
            if cc > 1       %if current cell is a duaghter
                hold on;
                plot(x(1),y(1), 'ko', 'markerfacecolor', 'k')   %plot a dot at frame of automatically found mitosis
                firstMitosisFrame = gp_rc(cc,1);
            end
            if selectmode | sample
                text(x(1),y(1)+0.2,num2str(cellid));
            end
        end
        annotation('textbox',[0.01,0.01,0.1,0.1],'string','','edgecolor','none');   %this and delete(findall()) required to avoid mouse sensitivity bug
        title(['Group ',num2str(tc_sel(counter))]);                
        if drugtime
            line([drugtime drugtime],[ymin ymax],'Color','r');  %add vertical line when drug is added
        end
        %%%% selectmode: view and choose cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 'enter': see next cell group
        % click on a trace: 'y','enter': saves cell id
        % 'i','enter': enter image viewing mode
        %   * must 'delete' previous datatips prior to entering image mode.
        %   - click on a trace to see image. arrow keys update next point.
        % 'd','enter': return to data viewing mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if selectmode
            dcm_obj = datacursormode(gcf);
            datacursormode on;
            response = 'on';
            while response              %get next graph: 'enter'
                set(dcm_obj,'Updatefcn',@datatip_traceid);
                response = input('','s');
                if strcmp(response,'y')                     %save traceid: 'y','enter'
                    figure(gcf);
                    nextindex = find(~savetracelist,1);
                    savetracelist(nextindex) = str2num(traceid);
                    ['saving cell ',traceid]
                end
                if strcmp(response,'i')                     %enter image mode: 'i','enter'
                    figure(gcf);
                    delete(findall(gcf,'Type','hggroup'));  %somehow prior annotations make plot mouse sensitive
                    set(dcm_obj,'Updatefcn',@datatip_image);
                    while ~strcmp(response,'d')             %exit image mode:'d','enter'
                        response = input('','s');
                    end
                    delete(findall(gcf,'Type','hggroup'));
                    close(gcf+1);       %close the image window
                end
                figure(gcf);
            end
        end
    end
end
%close(gcf);
if selectmode
    tracelist = savetracelist(find(savetracelist));
    save([savedir,moviename,'_goodcells_temp.mat'],'tracelist');
end
