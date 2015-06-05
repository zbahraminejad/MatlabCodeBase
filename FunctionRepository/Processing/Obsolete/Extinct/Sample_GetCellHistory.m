%% get red level
%Mingyu
%% preparation
% change cwd and tc
%clear;close all;clc
%cwd='F:\s1\s2\data\';
cwd='/Users/Mingyu/MATLAB/2012-01-15_notreat';
%this code has been updated to be more stringent in the cells it plots 

%load([cwd,'/data/5_12_0_alldata'],'bestsp','best_rc')
load([cwd,'/ring-2to0/5_12_1_alldata_5org6new7orgring8newring9intorg10intnew_raw'],'bestsp','best_rc')
%load([cwd,'/orgdata/5_12_0_alldata_5son6mctd7snn_raw'],'bestsp','best_rc')
load([cwd,'/data/cytodat'],'cytodat')
%numframes=133;  %(NOTE)
numframes=120; %Mingyu
framesperhr=5;


minratios=[];
%tracks=sort([47,12,129,90,98,42,163,88]);
%tracks=[5]
tracks=[55];
i=1;
while i<=size(tracks,2)
    if best_rc(tracks(i),2)~=best_rc(tracks(i),5)
        %[tracks(i),best_rc(tracks(i),5),best_rc(tracks(i),2)] %debug readout
        tracks=[tracks best_rc(tracks(i),2)];
    end
    i=i+1;
end
tracks=sort(tracks);
best_rc=best_rc(tracks,:);
for f=1:size(bestsp,3)
    ftracks=tracks(tracks<=size(bestsp{f},1));
    bestsp{f}=bestsp{f}(ftracks,:);
end

%% loops through each frame for all the cells

angieratio=-10000*ones(size(bestsp{end},1),size(bestsp,3));
cdt1=-10000*ones(size(bestsp{end},1),size(bestsp,3));

for f=1:size(angieratio,2)
    tempcell=find(bestsp{f}(:,1)~=0);
    cytoconc=(cytodat-bestsp{f}(tempcell,10))/
    angieratio(tempcell,f)=bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);
    %angieratio(tempcell,f)=bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);
    %angieratio(tempcell,f)=bestsp{f}(tempcell,5);  % calculating just nucleus
    %angieratio(tempcell,f)=bestsp{f}(tempcell,7);  % calculating just ring
    %angieratio(tempcell,f)=cytodat(f);  % calculating full cytosol
    %cdt1(tempcell,f)=bestsp{f}(tempcell,6);  %add data to angieratio. col 5 is reporter's intensity. col 3 is nuclear marker's intensity
end

%% remove negative values and assign them to be 0

angieratio(angieratio<0 & angieratio~=-10000)=0;
for track=1:size(bestsp{end},1)
    atrack=angieratio(track,:);
    atrack(isnan(atrack))=max(atrack);
    angieratio(track,:)=atrack;
    %angieratio(track,:)=atrack/max(atrack); % calculating only nuc or ring
end
%cdt1(cdt1<0 & cdt1~=-10000)=0

%% grouping
GP=zeros(size(best_rc,1),1);
counter=0;
for cc=1:size(GP,1)
    if best_rc(cc,2)~=best_rc(cc,5)
        
        %GP(cc)=GP(best_rc(cc,2));
        tracksinit=find(tracks==best_rc(cc,2),1);
        GP(cc)=GP(tracksinit);
        
        %fix angieratio
        angieratiomtx=angieratio([tracksinit,cc],:);
        %cdt1mtx=cdt1([tracksinit,cc],:);
        
        existlo=sum(angieratiomtx>-10000)==2;
        gpsplit=find(existlo&(diff(angieratiomtx,1,1)==0),1,'last');
        if (~isempty(gpsplit))&&(gpsplit<best_rc(cc,3))
            gpexist=find(existlo,1,'first');
            
            angieratio(tracksinit,gpexist:gpsplit)=mean(angieratiomtx(:,gpexist:gpsplit));
            angieratio(cc,1:gpsplit)=-10000;
            %cdt1(tracksinit,gpexist:gpsplit)=mean(cdt1mtx(:,gpexist:gpsplit));
            %cdt1(cc,1:gpsplit)=-10000;
            
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


tc_sel=[1:max(GP)];
trajectories = zeros(length(tc_sel),numframes);
n_tc = 0;


j=0;
for tc=tc_sel%:max(GP)% tc=tc_sel  %tc counts each group
    gp_rc0=best_rc(GP==tc,:);
    gp_rc0(:,6)=mod(1:size(gp_rc0,1),7)+1;  %to loop through the 7 colors and give each cell a color in the 6th col of gp_rc0
    gp_rc=gp_rc0;
    
    %%  make selected cells more stringent
    finalframecells=gp_rc0(:, 3)==numframes;  %only include cells that exist at final frame; the end frame for each cell is stored in col 3 of gp_rc
    gp_rc=gp_rc0(finalframecells, :);  %overwrite gp_rc to only include cells that exist at final frame
    %%%%
       

    
    %%%  make selected cells more stringent  (angie's sensor)
    %yfp=angieratio(gp_rc(:,5),:);
    
    
    %IF=cdt1(gp_rc(:,5),:);
    %intv=max([yfp(:); IF(:)]);
    %intv=max([yfp(:)]);
    
    
    %% only plot long tracks
    if (max(gp_rc(:,3))-min(gp_rc(:,1))+1)>=0  % *1 is the fraction of the length of the movie that at least one member of the group must exist for 
        
      
        
        j=j+1;
        %figure(ceil(j/1)); 
        set(gcf,'color','w')
        %subplot(1,1,mod(j-1,1)+1);  %for each subplot
        
        %set(gca,'YTick',0.5:0.5:2,'XTick', 0:framesperhr*4:numframes)
        
        %plot intensities:
        for cc=1:size(gp_rc,1)  %for each cell
            tempcolor=colors(gp_rc(cc,6));  %get color info from 6th col of gp_rc0
            
            angieratioeachcell=angieratio(find(tracks==gp_rc(cc,5),1),:);
            %cdt1eachcell=cdt1(gp_rc(cc,5),:);  %col 5 of gp_rc contains the # of each cell
            
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
            
            y3=angieratioeachcell(angieratioeachcell>-10000);        %y3 contains intensity thru time of fluorescent protein you wish to plot
            %y3=smooth(y3);
            
            %l3=line(x3(17:94),y3(17:94),'color','b','linewidth',4);
            
            l3=line(x3(1:end-1),y3(1:end-1),'color','r','linewidth',1);
            %l3=line(x3(1:end-1),y3(1:end-1),'color',tempcolor,'linewidth',4);  %this plots the YFP intensity through time 
            hold on
            %text(x3(1)+2,y3(1)+2,num2str(gp_rc(cc,5)),'horizontalalignment','center','fontsize',10)
            g=get(gca,'Position');      %[left fig edge to left plot edge, fig bottom to plot bottom, plot width, plot height]
            xl=get(gca,'XLim');
            yl=get(gca,'YLim');
            posx1=((x3(1)-xl(1))/range(xl))*g(3); % position relative to lower left corner of figure
            posy1=((y3(1)-yl(1))/range(yl))*g(4); % position relative to lower left corner of figure
            %annotation('textbox',[g(1)+posx1,g(2)+posy1,g(3)/10,g(4)/10],'string',num2str(gp_rc(cc,5)),'edgecolor','none');
            

          
 
            %% extract the first frame of mitosis and fit the slopes
            daughter_logic=best_rc(:,2)==gp_rc(cc,5);  %find identity of the daughter that is born from the current cell
            daughter_mitosis_frame=best_rc(daughter_logic,1);  %find the frame of mitosis for that daughter (there will only be one # if not mitosis); %the first element of the vector is 1 bc it's from itself   
            firstMitosisFrame=min(daughter_mitosis_frame(daughter_mitosis_frame~=1));
            mingyumitosis=mingyuindex(firstMitosisFrame);
            if mingyumitosis==0
                mingyumitosis=1;
            end
            plot(firstMitosisFrame, y3(mingyumitosis), 'ko', 'markerfacecolor', 'k')  %plot a dot at frame of automatically found mitosis
 
        end
        
        %xlim([17 94]);
        %title('med(nucleus)','fontsize',14);
        %xlabel('Frame (every 12 minutes)');
        %ylabel('Signal (YFP)');
        llim=0.12;
        ulim=0.36;
        %llim=0.2;
        %ulim=0.53;
        %range=ulim-llim;
        %ylim([llim ulim]);
        %set(gca,'yTick',llim:range/5:ulim)
        
    end
end

[cellname,fr]=find (cdt1>.9);  %find the index of cells with brightest EdU (>0.9)
%ylabel('Log (EdU signal:background)')

%save 3_2_0_non-mitosing_trajectories trajectories %