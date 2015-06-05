function sample_cell()
clear all;close all; clc;
cd('/Users/Steve/Dropbox/Matlab Files/From Mingyu/Functions');
path = '/Volumes/CAPPELL-2TB/20130323-MCF10A-EGF/';
savedir = '/Users/Steve/Documents/Meyer_Lab/Data/Fucci Sensors/Movies/';
datadir = ([path,'Data/']);
cpdir = [path,'CroppedProcessed/'];
%%% set average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sample specification %%%%%%%%%%%
row=4;
col=7;
site=1;
track=939;
signal_to_plot=2;  %% 1 for signal1, 2 for signal 2, or 3 for signal 1 and 2

%%% Setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,shot,'_alldata'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
windowsize = 20;
windowhalf = windowsize/2;

windowsize_final = 6;
windowhalf_final = windowsize_final/2;
%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[height,width]=size(tempread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startFrame=best_rc(track,1);
endFrame=best_rc(track,3);

SF=50;%best_rc(track,1);
EF=175;%best_rc(track,3);

frames=SF:EF;
movie_length=length(SF:EF);
%%% Get Signal traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal1 = -10000*ones(1,EF);
signal2 = -10000*ones(1,EF);
for ff=startFrame:endFrame
    signal1(ff) = bestsp{ff}(track,7)./bestsp{ff}(track,5);
    signal2(ff) = bestsp{ff}(track,6);
    signal2_area(ff) =bestsp{ff}(track,6).*bestsp{ff}(track,4);
end

%%% Test plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1),
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[800 scrsz(4)/1 scrsz(3)/2.3 scrsz(4)/3.5],'Color','w');
for kk=1:size(signal2,1);
    signal1_smooth(kk,:)=smooth(signal1(kk,:),1);
    signal2_smooth(kk,:)=smooth(signal2_area(kk,:),5);
end
background_to_add=0.01; %%For normalization of Geminin signal
signal2_norm=normalizeMyTracesGeminin_alt2(signal2_smooth,background_to_add);
signal2_norm_2=normalizeMyTracesGeminin_alt2(signal2_area,background_to_add);
APC=E3_activity_from_trace(signal2_norm);APC_last=APC(end);
signal2_norm=normalizeMyTracesGeminin_alt2(signal2_smooth,background_to_add);
% APC=[APC APC_last];
xaxis=0:0.2:size(signal2,2)/(5)-0.2; %make xaxis bases on time (hr)
subplot(1,2,1),hold on;
[haxes,hline1,hline2] = plotyy(xaxis,signal2_norm-background_to_add,xaxis,APC);
axes(haxes(1));axis([0 xaxis(end) 0 1]);
set(gca,'Box','off','YAxisLocation','left','YColor',[0.77 0 0],'YTick',0:0.2:1,'FontSize',16);
set(hline1,'color',[0.77 0 0],'linewidth',2);
set(get(haxes(1),'Ylabel'),'String','Geminin Levels (RFU)','FontSize',16)
    
axes(haxes(2));axis([0 xaxis(end) 0 1]);
set(gca,'YAxisLocation','right','YColor',[0 0 1],'YTick',0:0.2:1,'FontSize',16);
set(hline2,'color',[0 0 1],'linewidth',2);
set(get(haxes(2),'Ylabel'),'String','Relative APC activity','FontSize',16)
xlabel('Time (hr)','FontSize',16);

subplot(1,2,2),hold on;
[haxes,hline1,hline2] = plotyy(xaxis,signal1_smooth,xaxis,signal2_norm-background_to_add);
axes(haxes(1));axis([0 xaxis(end) 0.2 2]);
set(gca,'Box','off','YAxisLocation','left','YColor',[0 0.68 0],'YTick',0.2:0.2:2,'FontSize',16);
set(hline1,'color',[0 0.68 0],'linewidth',2);
set(get(haxes(1),'Ylabel'),'String','Relative CDK2 activity','FontSize',16)
    
axes(haxes(2));axis([0 xaxis(end) 0 1]);
set(gca,'YAxisLocation','right','YColor',[0.77 0 0],'YTick',0:0.2:1,'FontSize',16);
set(hline2,'color',[0.77 0 0],'linewidth',2);
set(get(haxes(2),'Ylabel'),'String','Geminin Levels (RFU)','FontSize',16)
xlabel('Time (hr)','FontSize',16);

%%% Make Movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aviobj=avifile([savedir,shot,'_',num2str(track),'_signal',num2str(signal_to_plot),'.avi'],'fps',7);
hh=figure(2);
set(gcf,'Position',[200 scrsz(4)/2.3 scrsz(3)/4.4 scrsz(4)/1.5],'Color','w');
arrow1=[0.57 0.53];arrow2=[0.79 0.74];arrow3=[0.31 0.31];arrow4=[0.43 0.40];arrow5=[0.48 0.48];
for f=SF:EF
    disp(f);
    %%% Read in pre-cropped images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    NEs_or=single(imread([cpdir,shot,'_nucedge&ring_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    CEs_or=single(imread([cpdir,shot,'_geminin_',num2str(f),'.tif']));
    
    %%%%%  Focus on track   %%%%%%%%%%%%%%
    cx=int16(bestsp{f}(track,2));
    cy=int16(bestsp{f}(track,1));
    miny=cy-(windowhalf-1)*nucr; maxy=cy+windowhalf*nucr;   
    minx=cx-(windowhalf-1)*nucr; maxx=cx+windowhalf*nucr;   
    %%%%%%%
    if minx<1
        minx=1; maxx=1+windowsize*nucr;
    end
    if miny<1
        miny=1; maxy=1+windowsize*nucr;
    end
    if maxx>height
        maxx=height; minx=height-windowsize*nucr;
    end
    if maxy>width
        maxy=width; miny=width-windowsize*nucr;
    end
    %%%%%%% Crop Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_large=DAs_or(minx:maxx,miny:maxy);     
    NEs_large=NEs_or(minx:maxx,miny:maxy);    
    REs_large=REs_or(minx:maxx,miny:maxy);     
    CEs_large=CEs_or(minx:maxx,miny:maxy);     
    
    REs_large=mat2gray(REs_large);      %normalize signal1 channel
    CEs_large=mat2gray(CEs_large);      %normalize signal2 channel
    NEs_large=mat2gray(NEs_large);      %normalize nucleus channel
    DAs_nucleus=mat2gray(DAs_large);  %normalize nucleus channel
    
    %%% Recrop
    crop_factor_1=(windowhalf-1)*nucr + windowhalf*nucr + 1;
    crop_factor_final=(windowhalf_final-1)*nucr + windowhalf_final*nucr + 1;
    start_row=((crop_factor_1-crop_factor_final)/2)+1;
    end_row=start_row+crop_factor_final-1;
    REs_large=REs_large(start_row:end_row,start_row:end_row);
    CEs_large=CEs_large(start_row:end_row,start_row:end_row);
    NEs_large=NEs_large(start_row:end_row,start_row:end_row);
    DAs_nucleus=DAs_nucleus(start_row:end_row,start_row:end_row);
    %%% Make small box with just Nucleus signal, surrounded by white line %
    size_of_white_bar=1;
    white_bar_nucleus_bottom=ones(size(DAs_nucleus,1),size_of_white_bar);
    white_bar_nucleus_bottom=single(white_bar_nucleus_bottom);
    DAs_nucleus_longer=[white_bar_nucleus_bottom DAs_nucleus];

    %%% Make large white box with small nucleus box in the bottom right corner
    nucleus_box_large=zeros(size(DAs_nucleus));nucleus_box_large=single(nucleus_box_large);
    nucleus_box_final=[nucleus_box_large, DAs_nucleus_longer];
    
    %%% Add the small nucleus box to the bottom right corner of signal 1 box    
    signal1_window=[REs_large, DAs_nucleus_longer];
    
    %%% Add the small nucleus box to the bottom right corner of signal 2 box 
    signal2_window=[CEs_large, DAs_nucleus_longer];
    %%% Add frame to movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if signal_to_plot==1;
        subplot(5,1,4:5);hold on
        plot(xaxis(startFrame:f-SF+1),signal1_smooth(SF:f),'color',[0 0.68 0],'linewidth',2);
        xlim([0 movie_length/5]);ylim([0.2 2]);
        set(gca,'YTick',0.2:0.2:2,'FontSize',16);
        ylabel('Relative CDK2 activity','FontSize',16);xlabel('Time (hr)','FontSize',16);
        %%%%%%
        tempframe=imadjust(nucleus_box_final);
        tempframe(:,:,2)=imadjust(signal1_window);
        tempframe(:,:,3)=imadjust(nucleus_box_final);
    
        figure(2),
        subplot(5,1,1:3);hold off
        imshow(imresize(tempframe,4));
        text(800,890,'H2B','Color','w','FontSize',16,'FontWeight','bold');
        text(25,890,'DHB-Ven','Color','w','FontSize',16,'FontWeight','bold');
        text2=annotation('arrow',arrow1,arrow2,'color','w','LineWidth',2);
%         aviobj=addframe(aviobj,getframe(gcf));
    end
    
    if signal_to_plot==2;
        subplot(5,1,4:5);hold on
        [haxes,hline1,hline2] = plotyy(xaxis(startFrame:f-SF+1),signal2_norm_2(SF:f)-background_to_add,xaxis(startFrame:f-SF+1),APC(SF:f));
        axes(haxes(1));axis([0 movie_length/5 0 1]);
        set(gca,'Box','off','YAxisLocation','left','YColor',[0.77 0 0],'YTick',0:0.2:1,'FontSize',16);
        set(hline1,'color',[0.77 0 0],'linewidth',2);
        set(get(haxes(1),'Ylabel'),'String','Geminin Levels (RFU)','FontSize',16)
    
        axes(haxes(2));axis([0 movie_length/5 0 1]);
        set(gca,'YAxisLocation','right','YColor','b','YTick',0:0.2:1,'FontSize',16);
        set(hline2,'color','b','linewidth',2);
        set(get(haxes(2),'Ylabel'),'String','Relative APC activity','FontSize',16)
        xlabel('Time (hr)','FontSize',16);
%         plot(xaxis(startFrame:f-SF+1),signal2_norm_2(SF:f),'color',[0.77 0 0],'linewidth',2);
%         xlim([0 movie_length/5]);ylim([0 1]);
%         set(gca,'YTick',0:0.2:1,'FontSize',16);
%         ylabel('Geminin Levels (RFU)','FontSize',16);xlabel('Time (hr)','FontSize',16);

        %%%%%%%
        tempframe=signal2_window;
        tempframe(:,:,2)=nucleus_box_final;
        tempframe(:,:,3)=nucleus_box_final;
    
        figure(2),
        subplot(5,1,1:3);hold off
        imshow(imresize(tempframe,4));
        text(350,-20,'H2B','Color','k','FontSize',16,'FontWeight','bold');
        text(60,-20,'mChy-Geminin','Color','k','FontSize',16,'FontWeight','bold');
        if f>=79
            annotation('arrow',arrow3,arrow4,'color','k','LineWidth',2);
            text(85,375,'Mitosis','Color','k','FontSize',16,'FontWeight','bold');
        end
        if f>=106
            annotation('arrow',arrow5,arrow4,'color','k','LineWidth',2);
            text(200,375,'G1/S','Color','k','FontSize',16,'FontWeight','bold');
        end
%         text2=annotation('arrow',arrow1,arrow2,'color','w','LineWidth',2);
        aviobj=addframe(aviobj,getframe(gcf));
    end
    
    if signal_to_plot==3;
        subplot(5,1,4:5);hold on
        [haxes,hline1,hline2] = plotyy(xaxis(startFrame:f-SF+1),signal1_smooth(SF:f),xaxis(startFrame:f-SF+1),signal2_norm(SF:f)-background_to_add);
        axes(haxes(1));axis([0 movie_length/5 0.2 2]);
        set(gca,'Box','off','YAxisLocation','left','YColor',[0 0.63 0],'YTick',0.2:0.2:2,'FontSize',16);
        set(hline1,'color',[0 0.63 0],'linewidth',2);
        set(get(haxes(1),'Ylabel'),'String','Relative CDK2 activity','FontSize',16)
    
        axes(haxes(2));axis([0 movie_length/5 0 1]);
        set(gca,'YAxisLocation','right','YColor',[0.77 0 0],'YTick',0:0.2:1,'FontSize',16);
        set(hline2,'color',[0.77 0 0],'linewidth',2);
        set(get(haxes(2),'Ylabel'),'String','Geminin Levels (RFU)','FontSize',16)
        xlabel('Time (hr)','FontSize',16);
        
        %%%%%%%
        tempframe=imadjust(signal2_window);
        tempframe(:,:,2)=imadjust(signal1_window);
        tempframe(:,:,3)=imadjust(nucleus_box_final);
    
        figure(2),
        subplot(5,1,1:3);hold off
        imshow(imresize(tempframe,4));
        text(800,890,'H2B','Color','w','FontSize',16,'FontWeight','bold');
        text(25,890,'mChy-Geminin','Color',[0.77 0 0],'FontSize',16,'FontWeight','bold');
        text(300,890,'DHB-Ven','Color',[0 0.68 0],'FontSize',16,'FontWeight','bold');
        text2=annotation('arrow',arrow1,arrow2,'color','w','LineWidth',2);
%         aviobj=addframe(aviobj,getframe(gcf));
    end
end
  aviobj=close(aviobj);
%   close(gcf)
end