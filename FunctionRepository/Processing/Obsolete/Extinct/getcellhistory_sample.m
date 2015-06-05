%% get red level
%% preparation
% change cwd and tc
clear;close all;clc
%cwd='F:\s1\s2\data\';
cwd='h:\Documents\Timescape\20120807_12drugs';
%this code has been updated to be more stringent in the cells it plots 
%%%%%%%%%%%% Experiment Settings %%%%%%%%
load([cwd,'\Data\4_12_1_alldata'],'bestsp','best_rc')
numframes=208; %Mingyu
framesperhr=5;
cdtplot = 0;
%%%%%%%%%%%% Set drug spike time %%%%%%%%%
SP = 90;
%%%%%%%%%%%% Set goups to plot %%%%%%%%%%%%
tc_sel = [68 66 37];
ymax = 2; ymin = 0.3;
minratios=[];

%% loops through each frame for all the cells

angieratio=-10000*ones(size(bestsp{end},1),size(bestsp,3));
cdt1=-10000*ones(size(bestsp{end},1),size(bestsp,3));

for f=1:size(angieratio,2)
    tempcell=find(bestsp{f}(:,1)~=0);
    angieratio(tempcell,f)=bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);  %change which fluorescent channel you plot here!  col 6 is antibody/EdU intensity; divide 7th col by 5th col to get cyto:nuc ratio of angie's sensor
    cdt1(tempcell,f)=bestsp{f}(tempcell,6);  %add data to angieratio. col 5 is reporter's intensity. col 3 is nuclear marker's intensity
end

%% remove negative values and assign them to be 0
angieratio(angieratio<0 & angieratio~=-10000)=0;
cdt1(cdt1<0 & cdt1~=-10000)=0;
%% grouping
GP=zeros(size(best_rc,1),1);
counter=0;
for cc=1:size(GP,1)
    if best_rc(cc,2)~=best_rc(cc,5)
        %split exits
        %group it
        GP(cc)=GP(best_rc(cc,2));
        %fix angieratio
        angieratiomtx=angieratio([best_rc(cc,2),cc],:);
        cdt1mtx=cdt1([best_rc(cc,2),cc],:);
        
        existlo=sum(angieratiomtx>-10000)==2;
        gpsplit=find(existlo&(diff(angieratiomtx,1,1)==0),1,'last');
        if (~isempty(gpsplit))&&(gpsplit<best_rc(cc,3))
            gpexist=find(existlo,1,'first');
            
            angieratio(best_rc(cc,2),gpexist:gpsplit)=mean(angieratiomtx(:,gpexist:gpsplit));
            angieratio(cc,1:gpsplit)=-10000;
            cdt1(best_rc(cc,2),gpexist:gpsplit)=mean(cdt1mtx(:,gpexist:gpsplit));
            cdt1(cc,1:gpsplit)=-10000;
            
            best_rc(cc,1)=(gpsplit+1);
        end
    else
        counter=counter+1;
        GP(cc)=counter;
    end
end

%% plotting
counter=0;
colors='ygbmkrc';
cmap=jet(128);cmap=[0,0,0;cmap];
differenceyfp=diff(cdt1,1,2);
trajectories = zeros(length(tc_sel),numframes);
n_tc = 0;
figure(1);
set(gcf,'color','w');
for counter=1:length(tc_sel)%:max(GP)% tc=tc_sel  %tc counts each group
    
    subaxis(1,length(tc_sel),counter, 'SH', 0.04, 'SV', 0.1,'Margin',0.05);
    tc = tc_sel(counter);
    gp_rc0=best_rc(GP==tc,:);
    gp_rc0(:,6)=mod(1:size(gp_rc0,1),7)+1; 
    gp_rc=gp_rc0;
    
    %%%%  make selected cells more stringent
    finalframecells=gp_rc0(:, 3)==numframes;
    gp_rc=gp_rc0(finalframecells, :); 
    %%%  make selected cells more stringent  (angie's sensor)
    yfp=angieratio(gp_rc(:,5),:);
    intv=max([yfp(:)]);
    
    
    %% only plot long tracks
        set(gca,'YTick',0.5:0.5:2,'XTick', 0:framesperhr*8:numframes)
        
        %plot intensities:
        for cc=1:size(gp_rc,1)  %for each cell      
            tempcolor=colors(gp_rc(cc,6));  %get color info from 6th col of gp_rc0
            angieratioeachcell=angieratio(gp_rc(cc,5),:);  %col 5 of gp_rc contains the # of each cell
            cdt1eachcell=cdt1(gp_rc(cc,5),:);  %col 5 of gp_rc contains the # of each cell
            
            x3=gp_rc(cc,1):gp_rc(cc,3);        %x3 is start frame through end frame for cell cc
       
            mingyuindex=angieratioeachcell>-10000;
            mingyuindex=single(mingyuindex);
            m=0;
            for a=1:size(mingyuindex,2)
                if mingyuindex(a)
                    m=m+1;
                    mingyuindex(a)=m;
                end
            end
            
            y3=angieratioeachcell(angieratioeachcell>-10000);
            y3=smooth(y3);
            l3=line(x3(1:end-1),y3(1:end-1),'color',tempcolor);
            if cdtplot
                y4=cdt1eachcell(cdt1eachcell>-10000);        %y4 contains intensity thru time 
                y4=mat2gray(y4)+0.5;                             %normalize cdt1 to be between 0 and 1, then add 0.5 so you can plot the Cdt1 on same y axis as angie's sensor
                l4=line(x3(1:end-1),y4(1:end-1),'color',tempcolor, 'linestyle', ':', 'linewidth', 1.2);  %this plots the CFP intensity through time 
            end
            g=get(gca,'Position');      %[left fig edge to left plot edge, fig bottom to plot bottom, plot width, plot height]
            xl=get(gca,'XLim');
            yl=get(gca,'YLim');
            posx1=((x3(1)-xl(1))/range(xl))*g(3); % position relative to lower left corner of figure
            posy1=((y3(1)-yl(1))/range(yl))*g(4); % position relative to lower left corner of figure
            annotation('textbox',[g(1)+posx1,g(2)+posy1,g(3)/10,g(4)/10],'string',num2str(gp_rc(cc,5)),'edgecolor','none');
            
            
            %% extract the first frame of mitosis and fit the slopes
            daughter_logic=best_rc(:,2)==gp_rc(cc,5);
            daughter_mitosis_frame=best_rc(daughter_logic,1);
            firstMitosisFrame=min(daughter_mitosis_frame(daughter_mitosis_frame~=1));
            mingyumitosis=mingyuindex(firstMitosisFrame);
            if mingyumitosis==0
                mingyumitosis=1;
            end
            hold on
            plot(firstMitosisFrame, y3(mingyumitosis), 'ko', 'markerfacecolor', 'k');
%             if cc==1
%                 ymax = max(y3); ymin = min(y3);
%             else
%                 if max(y3) > ymax
%                     ymax = max(y3);
%                 end
%                 if min(y3) < ymin
%                     ymin = min(y3);
%                 end
%             end
        end

        axis([0,numframes,ymin, ymax])
        %axis([0,size(bestsp,3),-0.1, intv*(size(gp_rc,1)/2+1+0.5)])
        axis xy
        
        %%%%%%%%%%%% Add vertical line and save image %%%%%%%%%%%%
        line([SP SP],[ymin ymax],'Color','r');
        title(['Group ',num2str(tc)]);
        hfig = figure(1);
        saveas(hfig,'h:\Downloads\TimeFig.jpg');
end

[cellname,fr]=find (cdt1>.9);  %find the index of cells with brightest EdU (>0.9)