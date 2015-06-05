%% preparation
clear;close all;clc;
path = 'h:\Documents\Timescape\20120807_12drugs\';
datadir = [path,'DataEval\'];
%%%%%%%%%%%% Experiment Settings %%%%%%%%
moviename = '3_1_1';
numframes = 6;
framesperhr = 5;
drugtime = 90;
%%%%%%%%%%%% Other Settings %%%%%%%%%
ymin = 0.3; ymax = 2;
historycutoff = 0.1;

load([datadir,moviename,'_alldata'],'bestsp','best_rc')
%% store signal values in angieratio and cdt1 matrices
angieratio=-10000*ones(size(bestsp{end},1),size(bestsp,3));
cdt1=-10000*ones(size(bestsp{end},1),size(bestsp,3));
for f=1:size(angieratio,2)
    tempcell=find(bestsp{f}(:,1)~=0); %ignore any cells that have negative coordinates (I'm not sure how negative coordinates would ever appear anyways...)
    angieratio(tempcell,f)=bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);  %cytosol/nuclear ratio for DHB (Angie's) sensor
    cdt1(tempcell,f)=bestsp{f}(tempcell,6);  %median nuclear Cdt1 value
end
% remove negative values and assign them to be 0
angieratio(angieratio<0 & angieratio~=-10000)=0;
cdt1(cdt1<0 & cdt1~=-10000)=0;

%% grouping
GP=zeros(size(best_rc,1),1); %no difference between size(bestsp,1) and size(best_rc,1)...
counter=0;
for cc=1:size(GP,1)
    if best_rc(cc,2)~=best_rc(cc,5) %split occurred (from mother to daughter)
        GP(cc)=GP(best_rc(cc,2));   %give the daughter the same GP value as mother (which occurred earlier, so already has a GP# assigned)
    else
        counter=counter+1;
        GP(cc)=counter;
    end
end

%% plotting
counter = 0;
colors = 'ygbmkrc';
tc_sel = [1:max(GP)];
trajectories = zeros(length(tc_sel),numframes);

j=0;
for tc = tc_sel%:max(GP)% tc=tc_sel  %tc counts each group
    gp_rc = best_rc(GP==tc,:); %store best_rc of all row(cell)#s in GP that have value 'tc'
    gp_rc(:,6) = mod(1:size(gp_rc,1),7)+1;  %to loop through the 7 colors and give each cell a color in the 6th col of gp_rc0
    if (max(gp_rc(:,3))-min(gp_rc(:,1))+1) >= (size(bestsp,3)*historycutoff)  % fraction of movie that at least one member of the group must exist for 
        j = j+1;
        figure(ceil(j/24)); 
        set(gcf,'color','w');
        subplot(4,6,mod(j-1,24)+1);  %for each subplot
        set(gca,'YTick',0.5:0.5:2,'XTick',0:framesperhr*4:numframes);

        %plot intensities:
        for cc = 1:size(gp_rc,1)  %for each cell      
            tempcolor = colors(gp_rc(cc,6));
            angieratioeachcell = angieratio(gp_rc(cc,5),:);  %get angieratio for each frame for this cell
            cdt1eachcell = cdt1(gp_rc(cc,5),:);  %same with cdt1
            x = gp_rc(cc,1):gp_rc(cc,3);        %x is start frame through end frame for cell cc
            
            mingyuindex = angieratioeachcell>-10000;
            mingyuindex = single(mingyuindex);
            m=0;
            for a = 1:size(mingyuindex,2)
                if mingyuindex(a)
                    m = m+1;
                    mingyuindex(a) = m;
                end
            end
            
            y = angieratioeachcell(angieratioeachcell>-10000);        %y contains intensity thru time of fluorescent protein you wish to plot
            line(x,y,'color',tempcolor);  %this plots the YFP intensity through time 
 
            %extract the first frame of mitosis and fit the slopes
            daughters = (best_rc(:,5)==gp_rc(cc,2));  %find all daughters that are born from the current cell
            daughter_mitosis_frame = best_rc(daughters,1);  %find the frame of mitosis for that daughter (there will only be one # if not mitosis); %the first element of the vector is 1 bc it's from itself   
            firstMitosisFrame = min(daughter_mitosis_frame(daughter_mitosis_frame~=1));
            mingyumitosis = mingyuindex(firstMitosisFrame);
            if mingyumitosis==0
                mingyumitosis = 1;
            end
            hold on
            plot(firstMitosisFrame, y(mingyumitosis), 'ko', 'markerfacecolor', 'k')  %plot a dot at frame of automatically found mitosis
        end

        axis([0,numframes,ymin,ymax]);
        axis xy
        if drugtime
            line([drugtime drugtime],[ymin ymax],'Color','r'); %add vertical line when drug is added
        end
        title(['Group ',num2str(tc)]);
    end
end