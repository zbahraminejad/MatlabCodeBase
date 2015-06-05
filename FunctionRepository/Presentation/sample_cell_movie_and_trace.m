function sample_cell_movie_and_trace()
%%%% Author: Steve Cappell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
%%% set average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sample specification %%%%%%%%%%%
row=1;
col=11;
site=1;
track=200;
signal_to_plot=1;  %% 1 for signal1, 2 for signal 2, or 3 for signal 1 and 2

%%% Setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,shot,'_alldata_HS_L05_595OfSegs_outer2'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
windowsize = 20;
windowhalf = windowsize/2;

windowsize_nucleus = 5;
windowhalf_nucleus = windowsize_nucleus/2;
%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[height,width]=size(tempread);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aviobj=avifile([datadir,shot,'_',num2str(track),'_signal',num2str(signal_to_plot),'.avi'],'fps',7);
startFrame=best_rc(track,1);
endFrame=best_rc(track,3);

SF=50;%best_rc(track,1);
EF=60;%best_rc(track,3);

frames=SF:EF;
movie_length=length(SF:EF);
%%% Get Signal traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal1 = -10000*ones(1,EF);
signal2 = -10000*ones(1,EF);
for ff=startFrame:endFrame
    signal1(ff) = bestsp{ff}(track,7)./bestsp{ff}(track,5);
    signal2(ff) = bestsp{ff}(track,6);
end

%%% Test plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1),
for kk=1:size(signal2,1);
    signal1_smooth(kk,:)=smooth(signal1(kk,:),5);
    signal2_smooth(kk,:)=smooth(signal2(kk,:),5);
end
signal2_norm=mat2gray(signal2_smooth); %normalize geminin signal
xaxis=0:0.2:size(signal2,2)/(5)-0.2; %make xaxis bases on time (hr)
[haxes,hline1,hline2] = plotyy(xaxis,signal1_smooth,xaxis,signal2_norm);
axes(haxes(1));axis([0 xaxis(end) 0.2 2]);
set(gca,'Box','off','YAxisLocation','left','YColor',[0 0.67 0],'YTick',0.2:0.2:2,'FontSize',16);
set(hline1,'color',[0 0.67 0],'linewidth',2);
set(get(haxes(1),'Ylabel'),'String','Relative CDK2 activity','FontSize',16)
    
axes(haxes(2));axis([0 xaxis(end) 0 1]);
set(gca,'YAxisLocation','right','YColor',[0.95 0 0],'YTick',0:0.2:1.015,'FontSize',16);
set(hline2,'color',[0.95 0 0],'linewidth',2);
set(get(haxes(2),'Ylabel'),'String','Geminin Levels (RFU)','FontSize',16)
xlabel('Time (hr)','FontSize',16);

%%% Make Movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh=figure(2);
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[800 scrsz(4)/2.3 scrsz(3)/4.4 scrsz(4)/1.5],'Color','w');
arrow1=[0.57 0.53];arrow2=[0.79 0.74];
for f=SF:EF
    disp(f);
    %%% Read in pre-cropped images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    NEs_or=single(imread([cpdir,shot,'_nucedge&ring_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    CEs_or=single(imread([cpdir,shot,'_cdt1_',num2str(f),'.tif']));
    
    %%%%%  Focus on track   %%%%%%%%%%%%%%
    cx=int16(bestsp{f}(track,2));
    cy=int16(bestsp{f}(track,1));
    miny=cy-(windowhalf-1)*nucr; maxy=cy+windowhalf*nucr;   miny_nucleus=cy-(windowhalf_nucleus-1)*nucr; maxy_nucleus=cy+windowhalf_nucleus*nucr;
    minx=cx-(windowhalf-1)*nucr; maxx=cx+windowhalf*nucr;   minx_nucleus=cx-(windowhalf_nucleus-1)*nucr; maxx_nucleus=cx+windowhalf_nucleus*nucr;
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
    %%%%%%%
    if minx_nucleus<1
        minx_nucleus=1; maxx_nucleus=1+windowsize_nucleus*nucr;
    end
    if miny_nucleus<1
        miny_nucleus=1; maxy_nucleus=1+windowsize_nucleus*nucr;
    end
    if maxx_nucleus>height
        maxx_nucleus=height; minx_nucleus=height-windowsize_nucleus*nucr;
    end
    if maxy_nucleus>width
        maxy_nucleus=width; miny_nucleus=width-windowsize_nucleus*nucr;
    end
    %%%%%%% Crop Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_large=DAs_or(minx:maxx,miny:maxy);     DAs_nucleus=DAs_or(minx_nucleus:maxx_nucleus,miny_nucleus:maxy_nucleus);
    NEs_large=NEs_or(minx:maxx,miny:maxy);     NEs_nucleus=NEs_or(minx_nucleus:maxx_nucleus,miny_nucleus:maxy_nucleus);
    REs_large=REs_or(minx:maxx,miny:maxy);     REs_nucleus=REs_or(minx_nucleus:maxx_nucleus,miny_nucleus:maxy_nucleus);
    CEs_large=CEs_or(minx:maxx,miny:maxy);     CEs_nucleus=CEs_or(minx_nucleus:maxx_nucleus,miny_nucleus:maxy_nucleus);
    
    REs_large=mat2gray(REs_large);      %normalize signal1 channel
    CEs_large=mat2gray(CEs_large);      %normalize signal2 channel
    DAs_nucleus=mat2gray(DAs_nucleus);  %normalize nucleus channel
    
    %%% Make small box with just Nucleus signal, surrounded by white line %
    size_of_white_bar=1;
    white_bar_nucleus_top=ones(size_of_white_bar,size(DAs_nucleus,2)+size_of_white_bar);
    white_bar_nucleus_top=single(white_bar_nucleus_top);
    white_bar_nucleus_bottom=ones(size(DAs_nucleus,1),size_of_white_bar);
    white_bar_nucleus_bottom=single(white_bar_nucleus_bottom);
    DAs_nucleus=[white_bar_nucleus_bottom DAs_nucleus];
    DAs_nucleus=[white_bar_nucleus_top;DAs_nucleus];
    
    %%% Make large white box with small nucleus box in the bottom right corner
    nucleus_box_large=zeros(size(DAs_large));nucleus_box_large=single(nucleus_box_large);
    nucleus_box_large_top=nucleus_box_large(1:size(DAs_large,1)-size(DAs_nucleus,1),1:size(DAs_large,2));
    nucleus_box_large_bottom=nucleus_box_large(1:size(DAs_nucleus,1),1:size(DAs_large,2)-size(DAs_nucleus,2));
    nucleus_box_large_bottom=[nucleus_box_large_bottom, DAs_nucleus];
    nucleus_box_final=[nucleus_box_large_top; nucleus_box_large_bottom];
    
    %%% Add the small nucleus box to the bottom right corner of signal 1 box
    REs_large_top=REs_large(1:size(REs_large,1)-size(DAs_nucleus,1),1:size(REs_large,2));
    REs_large_bottom=REs_large(size(REs_large,1)-size(DAs_nucleus,1)+1:size(REs_large,1),1:size(REs_large,2)-size(DAs_nucleus,2));
    REs_large_bottom=[REs_large_bottom, DAs_nucleus];
    signal1_window=[REs_large_top; REs_large_bottom];
    
    %%% Add the small nucleus box to the bottom right corner of signal 2 box 
    CEs_large_top=CEs_large(1:size(CEs_large,1)-size(DAs_nucleus,1),1:size(CEs_large,2));
    CEs_large_bottom=CEs_large(size(CEs_large,1)-size(DAs_nucleus,1)+1:size(CEs_large,1),1:size(CEs_large,2)-size(DAs_nucleus,2));
    CEs_large_bottom=[CEs_large_bottom, DAs_nucleus];
    signal2_window=[CEs_large_top; CEs_large_bottom];

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
        aviobj=addframe(aviobj,getframe(gcf));
    end
    
    if signal_to_plot==2;
        subplot(5,1,4:5);hold on
        plot(xaxis(startFrame:f-SF+1),signal2_norm(SF:f),'color',[0.77 0 0],'linewidth',2);
        xlim([0 movie_length/5]);ylim([0 1]);
        set(gca,'YTick',0:0.2:1,'FontSize',16);
        ylabel('Geminin Levels (RFU)','FontSize',16);xlabel('Time (hr)','FontSize',16);

        %%%%%%%
        tempframe=imadjust(signal2_window);
        tempframe(:,:,2)=imadjust(nucleus_box_final);
        tempframe(:,:,3)=imadjust(nucleus_box_final);
    
        figure(2),
        subplot(5,1,1:3);hold off
        imshow(imresize(tempframe,4));
        text(800,890,'H2B','Color','w','FontSize',16,'FontWeight','bold');
        text(25,890,'mChy-Geminin','Color','w','FontSize',16,'FontWeight','bold');
        text2=annotation('arrow',arrow1,arrow2,'color','w','LineWidth',2);
        aviobj=addframe(aviobj,getframe(gcf));
    end
    
    if signal_to_plot==3;
        subplot(5,1,4:5);hold on
        [haxes,hline1,hline2] = plotyy(xaxis(startFrame:f-SF+1),signal1_smooth(SF:f),xaxis(startFrame:f-SF+1),signal2_norm(SF:f));
        axes(haxes(1));axis([0 movie_length/5 0.2 2]);
        set(gca,'Box','off','YAxisLocation','left','YColor',[0 0.63 0],'YTick',0.2:0.2:2,'FontSize',16);
        set(hline1,'color',[0 0.63 0],'linewidth',2);
        set(get(haxes(1),'Ylabel'),'String','Relative CDK2 activity','FontSize',16)
    
        axes(haxes(2));axis([0 movie_length/5 0 1]);
        set(gca,'YAxisLocation','right','YColor',[0.77 0 0],'YTick',0:0.2:1.015,'FontSize',16);
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
        aviobj=addframe(aviobj,getframe(gcf));
    end
end
  aviobj=close(aviobj);
%   close(gcf)
end