%% preparation
% change cwd and tc
clear;close all;clc
%cwd='F:\s1\s2\data\';
path='h:\Documents\Timescape\20120807_12drugs\';
datadir = ([path,'DataEval\']);
%this code has been updated to be more stringent in the cells it plots
%%
row=[3];
for col=[1]
    for site=[1]
        movieName=[num2str(row),'_', num2str(col), '_', num2str(site)];  %[CHANGE]
        load([datadir, movieName, '_alldata'],'bestsp','best_rc')
        numframes=51;
        framesperhr=5;
        mitogenremovalframe=120;
        %% loops through each frame for all the cells
        angieratio=-10000*ones(size(bestsp{end},1),size(bestsp,3));
        cdt1=-10000*ones(size(bestsp{end},1),size(bestsp,3));
        for f=1:size(angieratio,2)
            tempcell=find(bestsp{f}(:,1)~=0);
            angieratio(tempcell,f)=bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);  %change which fluorescent channel you plot here!  col 6 is antibody/EdU intensity; divide 7th col by 5th col to get cyto:nuc ratio of angie's sensor
            %angieratio(tempcell,f)=bestsp{f}(tempcell,5);
            cdt1(tempcell,f)=bestsp{f}(tempcell,6);  %add data to angieratio. col 6 is Ab stain in far red.  called "Cdt1" for historical reasons
        end
        %% remove negative values and assign them to be 0
        angieratio(angieratio<0 & angieratio~=-10000)=0;
        cdt1(cdt1<0 & cdt1~=-10000)=0;
        
        %% loop through each cell and smoothen its trajectory
        
        % for c=1:size(angieratio,1)
        %     thiscell=angieratio(c,:);
        %     keep=thiscell~=-10000;
        %     smoothenedcell=smooth(thiscell(keep),5);  %default is to smoothen by averaging 5 time points
        %     angieratio(c,keep)=smoothenedcell;
        % end
        %
        % for c=1:size(cdt1,1)
        %     thiscell=cdt1(c,:);
        %     keep=thiscell~=-10000;
        %     smoothenedcell=smooth(thiscell(keep),5);  %default is to smoothen by averaging 5 time points; 1 is to not smoothen
        %     cdt1(c,keep)=smoothenedcell;
        % end
        
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
        colors='kgbmkrc';
        
        cmap=jet(128);cmap=[0,0,0;cmap];
        
        % %to test how bumpy the trajectories are
        % differenceyfp=diff(angieratio,1,2);
        % figure(21)
        % hist(differenceyfp(:),-60:.1:60)
        % xlim([-10 10])
        
        
        %to test how bumpy the trajectories are
        % differenceyfp=diff(cdt1,1,2);
        % figure(21)
        % hist(differenceyfp(:), -5000:10:5000)
        
        
        tc_sel=[1:max(GP)];
        trajectories = zeros(length(tc_sel),numframes);
        n_tc = 0;
        
        %stuff to store
        angieThruTime=zeros(numframes, 1000);
        cdt1ThruTime=zeros(numframes, 1000);
        frameOfMitosis=zeros(1, 1000);
        frameOfFirstMitosisAfterMitogenRemoval=zeros(1, 1000);
        groupNumber=zeros(1, 1000);
        %Abstain=zeros(1, 1000);
        
        
        j=0;
        for tc=1:max(GP)% tc=tc_sel  %tc counts each group
            gp_rc0=best_rc(GP==tc,:);
            gp_rc0(:,6)=mod(1:size(gp_rc0,1),7)+1;  %to loop through the 7 colors and give each cell a color in the 6th col of gp_rc0
            gp_rc=gp_rc0;
            
            %%  make selected cells more stringent
            %     finalframecells=gp_rc0(:, 3)==numframes;  %only include cells that exist at final frame; the end frame for each cell is stored in col 3 of gp_rc
            %     gp_rc=gp_rc0(finalframecells, :);  %overwrite gp_rc to only include cells that exist at final frame
            %     %%%%
            
            
            
            %     %%%  make selected cells more stringent  (angie's sensor)
            yfp=angieratio(gp_rc(:,5),:);
            deltayfp=diff(yfp,1,2);  %change in yfp intensity for each time step;
            deltayfp_bumpycells=(abs(deltayfp)>1.0)&(yfp(:,1:end-1)~=-10000)&(yfp(:,2:end)~=-10000); %if jump in yfp intensity >500, remove the cell. also account for the fact that the matrix was initialized at 10,000 so when a cell is born it will have a big jump
            smoothcells_logic=sum(deltayfp_bumpycells,2)==0;  %sum the rows and if it's equal to 0, then it's a smooth cell
            
            %     %%%
            
            %             %     %%%  make selected cells more stringent (Cdt1)
            %                 yfp=cdt1(gp_rc(:,5),:);
            %                 deltayfp=diff(yfp,1,2);  %change in yfp intensity for each time step;
            %                 deltayfp_bumpycells=(abs(deltayfp)>500) & (yfp(:,1:end-1)>exp(4))&(yfp(:,2:end)>exp(4)); %if jump in yfp intensity >300, remove the cell. also remove dim cdt1 cells (where signal is <exp(4)). also account for the fact that the matrix was initialized at 10,000 so when a cell is born it will have a big jump
            %                 smoothcells_logic=sum(deltayfp_bumpycells,2)==0;  %sum the rows and if it's equal to 0, then it's a smooth cell
            %                 gp_rc=gp_rc(smoothcells_logic,:);  %overwrite gp_rc with good cells
            %                 %%%
            
            
            
            
            
            %% only plot long tracks
            if (max(gp_rc(:,3))-min(gp_rc(:,1))+1)>=(size(bestsp,3)*1)  % *1 is the fraction of the length of the movie that at least one member of the group must exist for
                
                
                
                
                
                
                for cc=1:size(gp_rc,1)  %for each cell
                    disp(cc);
                    tempcolor=colors(gp_rc(cc,6));  %get color info from 6th col of gp_rc0
                    
                    
                    
                    angieratioeachcell=angieratio(gp_rc(cc,5),:);  %col 5 of gp_rc contains the # of each cell
                    cdt1eachcell=cdt1(gp_rc(cc,5),:);  %col 5 of gp_rc contains the # of each cell
                    
                    if gp_rc(cc,2) ~= gp_rc(cc,5)
                        missingangie=angieratio(gp_rc(cc,2), gp_rc(gp_rc(:,5)==gp_rc(cc,2),1):(gp_rc(cc,1)-1));
                        missingcdt1=cdt1(gp_rc(cc,2), gp_rc(gp_rc(:,5)==gp_rc(cc,2),1):(gp_rc(cc,1)-1));
%                        missingnuc=nucarea(gp_rc(cc,2), gp_rc(gp_rc(:,5)==gp_rc(cc,2),1):(gp_rc(cc,1)-1));
                        
                        angieratioeachcell(gp_rc(gp_rc(:,5)==gp_rc(cc,2),1):(gp_rc(cc,1)-1))=missingangie;
                        angieratio(gp_rc(cc,5),:)=angieratioeachcell;
                        
                        cdt1eachcell(gp_rc(gp_rc(:,5)==gp_rc(cc,2),1):(gp_rc(cc,1)-1))=missingcdt1;
                        cdt1(gp_rc(cc,5),:)=cdt1eachcell;
                        
%                         nucareaeachcell(gp_rc(gp_rc(:,5)==gp_rc(cc,2),1):(gp_rc(cc,1)-1))=missingnuc;
%                         nucarea(gp_rc(cc,5),:)=nucareaeachcell;
                    end
                    
                    
                    
                    
                    
                    %% extract the first frame of mitosis, *for both sisters*
                    
                    if gp_rc(cc,2)==gp_rc(cc,5)  %if the cell is the mother
                        
                        daughter_logic=best_rc(:,2)==gp_rc(cc,5);  %find identity of the daughter that is born from the current cell
                        daughter_mitosis_frame=best_rc(daughter_logic,1);  %find the frame of mitosis for that daughter (there will only be one # if not mitosis); %the first element of the vector is 1 bc it's from itself
                        firstMitosisFrame=min(daughter_mitosis_frame(daughter_mitosis_frame~=1));
                        
                    else  %if the cell is the daughter
                        firstMitosisFrame=gp_rc(cc,1);
                        gp_rc(cc,1)=gp_rc(gp_rc(:,5)==gp_rc(cc,2),1);  %assign first frame of daughter to be from the first frame of mother
                    end
                    
                    
                    
                    %% intereact with the figure to get frame of 1st S phase after black dot
                    %
                    %if trace is bad,  just hit enter.  if trace is good, hit 'y' and then enter
                    
                    %figure(tc)
                    %             subplot(2,1,2)
                    
                    %if ~isempty(firstMitosisFrame) %if a frame of mitosis is found
                    if ((firstMitosisFrame< (numframes)) + (firstMitosisFrame >gp_rc(cc,1)) + (firstMitosisFrame <gp_rc(cc,3))==3)  %if cell had mitosis 10 frames before end of movie
                        if gp_rc(cc,1)==1 && gp_rc(cc,3)==numframes  % if track is full length of movie
                            if smoothcells_logic(cc)==1  %if the cell is smooth, then plot it
                                
                                figure(tc), hold on
                                subplot(2,1,1), hold on
                                x3=gp_rc(cc,1):gp_rc(cc,3);        %x3 is start frame through end frame for cell cc
                                y3=angieratioeachcell(gp_rc(cc,1):gp_rc(cc,3));        %y3 contains intensity thru time of fluorescent protein you wish to plot
                                %y3=angieratioeachcell(angieratioeachcell>-10000);        %y3 contains intensity thru time of fluorescent protein you wish to plot
                                
                                %                 y3a=smooth(y3(1:firstMitosisFrame-gp_rc(cc,1)) );
                                %                 y3b=smooth(y3(firstMitosisFrame-gp_rc(cc,1)+1:end));
                                %                 y3=[y3a; y3b];
                                y3=smooth(y3);
                                
                                l3=line(x3(1:end-1),(y3(1:end-1)),'color',tempcolor);  %this plots the YFP intensity through time
                                %plot(firstMitosisFrame, y3(firstMitosisFrame), 'ko', 'markerfacecolor', 'k')  %plot a dot at frame of automatically found mitosis
                                ylim([.2 2.4])
                                title(['Group ',num2str(tc)])
                                plot([mitogenremovalframe mitogenremovalframe], [.25 2.25], 'k-')
                                
                                
                                
                                subplot(2,1,2), hold on
                                x3=gp_rc(cc,1):gp_rc(cc,3);        %x3 is start frame through end frame for cell cc
                                y4=cdt1eachcell(gp_rc(cc,1):gp_rc(cc,3));        %y3 contains intensity thru time of fluorescent protein you wish to plot
                                %y3=angieratioeachcell(angieratioeachcell>-10000);        %y3 contains intensity thru time of fluorescent protein you wish to plot
                                
                                %                 y3a=smooth(y3(1:firstMitosisFrame-gp_rc(cc,1)) );
                                %                 y3b=smooth(y3(firstMitosisFrame-gp_rc(cc,1)+1:end));
                                %                 y3=[y3a; y3b];
                                y4=smooth(y4);
                                
                                l4=line(x3(1:end-1),(y4(1:end-1)),'color',tempcolor);  %this plots the YFP intensity through time
                                plot([mitogenremovalframe mitogenremovalframe], [0 max(y4)], 'k-')
                                %plot(firstMitosisFrame, y4(firstMitosisFrame), 'ko', 'markerfacecolor', 'k')  %plot a dot at frame of automatically found mitosis
                                %ylim([.2 2.2])
                                
                                
                                %% interact with the plots
                                if ~isempty(firstMitosisFrame) %if a frame of mitosis is found
                                    [xcoord, ycoord]=getpts(gcf);
                                    xcoord=round(xcoord);
                                    
                                    if ~isempty(xcoord)  %if xcoord is not empty
                                        j=j+1;
                                        
                                        if xcoord(1)>numframes
                                            xcoord=1000;  %if cell doesn't enter s phase after mitosis, click in gray area to right of the plot; then change whatever # that is to 1000
                                        end
                                        
                                        angieThruTime(:,j)=y3';
                                        cdt1ThruTime(:,j)=y4';
                                        frameOfMitosis(j)=firstMitosisFrame;
                                        frameOfFirstMitosisAfterMitogenRemoval(j)=xcoord(1);
                                        groupNumber(j)=tc;
                                    end
                                end
                                close(tc)
                                
                                
                                

                            end
                        end
                        
                    end
                    
                    
                    
                end
                
                
                
                
            end
        end
        
        %remove extra zeros
        angieThruTime=angieThruTime(:,1:j);
        cdt1ThruTime=cdt1ThruTime(:,1:j);
        frameOfMitosis=frameOfMitosis(:,1:j);
        frameOfFirstMitosisAfterMitogenRemoval=frameOfFirstMitosisAfterMitogenRemoval(:,1:j);
        groupNumber=groupNumber(:,1:j);
        
        %        Abstain_trackedcells=Abstain(:,1:j);
        
        %        cdt1_logic=cdt1(:,end)>-10000;
        %        Abstain_allcells=cdt1(cdt1_logic, end);% store Ab intensities of all cells in final movie frame where intensity is not -10,000
        
        
        
        %save.  'cdt1' is the Abstain intensity for all cells in the field;
        %'Abstain' is the intensity for just the correctly tracked cells that have a mitosis
        
        save([datadir, movieName, '_dataForManuallySelectedCells_frameOfFirstMitosisAfterMitogenRemoval'],'angieThruTime','cdt1ThruTime', 'frameOfMitosis', 'frameOfFirstMitosisAfterMitogenRemoval', 'groupNumber')
        
        
        %save([cwd,'data\', movieName, '_dataForManuallySelectedCells'],'angieThruTime','frameOfMitosis', 'groupNumber')
        %save([cwd,'data\', movieName, '_dataForManuallySelectedCells'],'angieThruTime','frameOfMitosis', 'groupNumber', 'Abstain_trackedcells', 'Abstain_allcells')
        
        
        
    end
end


