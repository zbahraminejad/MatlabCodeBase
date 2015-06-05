row = 5; col = 2; site = 1;

% projectpath='D:\Documents\Projects\';
imagepath='G:\Michael\';
%imagepath='E:\';
experimentpath='20150401-CC-Diff\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
wellname = nameandsite(shot);
%shot=wellnum2str(row,col,site);
datadir=([imagepath,experimentpath,'Data\']);
separatedirectories=0;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    %rawdir=[imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
%     rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
    rawdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end
motheroption = 0;
daughteroption = 0;
IFoption = 0;
datafile =[datadir,'tracedata_',shot,'_nolink','.mat'];
load(datafile,'tracedata','genealogy','jitters');
[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother]=gathertracedata_mz_1(datadir,datafile,shot,motheroption,daughteroption,IFoption);

%%%%%% Set subtightplot settings %%%%%%%%%%%%%%%%%%%%%%%
gapcoordinates = [0.06 0.04];
heightmargin = 0.07;
widthmargin = 0.025;

%%%%%Get x and y coordinates%%%%%
traceid=[2276]; %180 228 350 %input cell id
fig = figure(1);
set(fig,'color','white','units','normalized','outerposition',[0 0 1 1]);
% subtightplot(2,3,[5,6],gapcoordinates,heightmargin,widthmargin)

time = 1:653;
xtime = time./5.5;
nucarea = tracedata(traceid,:,3);
pparg = tracedata(traceid,:,10);
ysig1 = smoothignorenans(pparg.*nucarea,6);
gem = tracedata(traceid,:,11);
ysig2 = smoothignorenans(gem,3);
CDK2nuc = tracedata(traceid,:,6);
CDK2cyto = tracedata(traceid,:,7);
ysig3 = smoothignorenans(CDK2cyto./CDK2nuc,6);
%%%%normalization
ysig2 = ysig2./max(ysig2);


ymin1 = 0; ymax1 = 150000; ystep1 = 50000;
ymin2 = 0; ymax2 = 3; ystep2 = 1;
tracewidth = 4;
[haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
axes(haxes(1));

axis([xtime(1) xtime(end) ymin1 ymax1]);
set(gca,'Box','off','YAxisLocation','left','YColor','k','YTick',ymin1:ystep1:ymax1);
set(hline1,'DisplayName',num2str(i),'color',[0 0.8 0],'linewidth',tracewidth); % was blue before 'b'
%gold = [0.9 .75 0]
axes(haxes(2));
axis([xtime(1) xtime(end) ymin2 ymax2]);
set(gca,'YAxisLocation','right','YColor','r','YTick',ymin2:ystep2:ymax2); % Ytick spacing default is --> ymin2:ystep2:ymax2
set(hline2,'DisplayName',num2str(i),'Color',[0.8 0 0],'linewidth',tracewidth);

trace3color = [0 0 0.8];
line(xtime,ysig3,'color',trace3color,'DisplayName',num2str(i),'linewidth',tracewidth);


hold on

% axes(haxes(1));

line([113/6 113/6] ,[0 max(ysig1)],'Color','k','linewidth',3,'linestyle','-');
line([385/6 385/6],[0 max(ysig1)],'Color','k','linewidth',3,'linestyle','--');
hold off
set(haxes(1),'FontName','Arial','FontSize',20);
xlabel(haxes(1),'Time (hours)','FontName','Arial','FontSize',26,'Units','Normalized','Position',[.5 -0.07]);
set(haxes(2),'FontName','Arial','FontSize',20);
axes(haxes(1))
dot1 = line('XData',xtime(1), 'YData',ysig1(1), 'Color','g','Marker','+', 'MarkerSize',14, 'LineWidth',10);
axes(haxes(2))
dot2 = line('XData',xtime(1), 'YData',ysig2(1), 'Color','r','Marker','+', 'MarkerSize',14, 'LineWidth',10);
dot3 = line('XData',xtime(1), 'YData',ysig3(1), 'Color','b','Marker','+', 'MarkerSize',14, 'LineWidth',10);

title('Time Traces','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
xlabel('Time (hours)','FontName','Arial','FontSize',26,'Units','Normalized','Position',[.5 -0.07]);
% ylabel('RFU','FontName','Arial','FontSize',26,'Units','Normalized','Position',[-.09 .5]);
% text(3, 25, '+Rosi', 'Color',[0 0 0],'FontSize',30, 'HorizontalAlignment','left', 'VerticalAlignment','top');
% text(48, 40, '+Insulin', 'Color',[0 0 0], 'FontSize',30, 'HorizontalAlignment','left', 'VerticalAlignment','top');

videofile ='C:\Users\teruellab\Desktop\example_4color_cell.avi';
% M=VideoWriter(videofile,'Uncompressed AVI');
M=VideoWriter(videofile,'Motion JPEG AVI');
M.FrameRate=26;
open(M);

 %setup to get max and min intensity 
channel = {'Cy5_','CFP_','YFP_','mCherry_'};
realdumb1=double(imread([rawdir,channel{2},num2str(200),'.tif']));
realdumb2=double(imread([rawdir,channel{3},num2str(200),'.tif']));
realdumb3=double(imread([rawdir,channel{2},num2str(200),'.tif']));

maskdumb=double(imread([maskdir,'nucedge_',num2str(200),'.tif']));
nucdumb = double(imread([rawdir,channel{1},num2str(200),'.tif']));
minInt1 = min(pparg);
maxInt1 = max(pparg)*0.75;
minInt3 = min(CDK2cyto);
maxInt3 = max(CDK2cyto)*0.6;
minInt2 = min(gem);
maxInt2 = max(gem)*0.75;

for frame = 1:653;  

allx=tracedata(:,frame,1)-jitters(frame,1); ally=tracedata(:,frame,2)-jitters(frame,2);
cx=int16(allx(traceid));
cy=int16(ally(traceid));

channel = {'Cy5_','YFP_','mCherry_','CFP_'};
bgperctile=10;
maxval=max(tracedata(traceid,:,6),[],2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real1=double(imread([rawdir,channel{2},num2str(frame),'.tif']));
real2=double(imread([rawdir,channel{3},num2str(frame),'.tif']));
real3=double(imread([rawdir,channel{4},num2str(frame),'.tif']));
mask=double(imread([maskdir,'nucedge_',num2str(frame),'.tif']));
nuc = double(imread([rawdir,channel{1},num2str(frame),'.tif']));

nucr = 8;
winrad=3*nucr;    %window radius: default 5
imagescale=4;  
%%%%% crop image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(real1);
%minpix=median(prctile(raw,0));  %sets absolute intensities. default 1
%maxpix=median(prctile(raw,100)); %default 99
minx=max([cx-winrad 1]);
maxx=min([cx+winrad width]);
miny=max([cy-winrad 1]);
maxy=min([cy+winrad height]);
cropmask=mask(miny:maxy,minx:maxx);
cropreal1=real1(miny:maxy,minx:maxx);
cropreal2=real2(miny:maxy,minx:maxx);
cropreal3=real3(miny:maxy,minx:maxx);
cropnuc = nuc(miny:maxy,minx:maxx);


cropmask1=imresize(cropmask,imagescale);
cropmask2=imresize(cropmask,imagescale);
cropmask3=imresize(cropmask,imagescale);

cropmask1(:,:,2)=mat2gray(imresize(cropreal1,imagescale),[minInt1 maxInt1]);
cropmask1(:,:,3)=0;
cropmask1(:,:,1)= imresize(cropmask,imagescale);

cropmask2(:,:,2)=0;
cropmask2(:,:,3)=imresize(cropmask,imagescale);
cropmask2(:,:,1)= mat2gray(imresize(cropreal2,imagescale),[minInt2 maxInt2]);

cropmask3(:,:,2)=0;
cropmask3(:,:,3)= mat2gray(imresize(cropreal3,imagescale),[minInt3 maxInt3]);
cropmask3(:,:,1)= imresize(cropmask,imagescale);

cropnuc =  mat2gray(imresize(cropnuc,imagescale));
tempnuc = cropnuc;

cropnuc(:,:,1) = tempnuc;
cropnuc(:,:,2) = tempnuc;
cropnuc(:,:,3) = tempnuc;


figure(1), subtightplot(2,3,1,gapcoordinates,heightmargin,widthmargin)
image(cropmask1);
set(gca,'YTick',[])
set(gca,'XTick',[])
title('PPARg','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
% axis square

figure(1), subtightplot(2,3,3,gapcoordinates,heightmargin,widthmargin)
image(cropmask3);
set(gca,'YTick',[])
set(gca,'XTick',[])
title('CDK2','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);

figure(1), subtightplot(2,3,2,gapcoordinates,heightmargin,widthmargin)
image(cropmask2);
set(gca,'YTick',[])
set(gca,'XTick',[])
title('Geminin','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);

figure(1) ,subtightplot(2,3,4,gapcoordinates,heightmargin,widthmargin)
image(cropnuc)
set(gca,'YTick',[])
set(gca,'XTick',[])
title('H2B','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
% axis square



figure(1)%,hold on, subtightplot(2,3,[5,6],gapcoordinates,heightmargin,widthmargin)
hold on
set(dot1,'XData',xtime(frame),'Ydata',ysig1(frame));
set(dot2,'XData',xtime(frame),'Ydata',ysig2(frame));
set(dot3,'XData',xtime(frame),'Ydata',ysig3(frame));
hold off
drawnow 
% axis square

f = getframe(fig);
writeVideo(M,f);
% addframe(videofile,fig.Number);

end
close(M);


%%
traceid = 180;
delay = 0.05;
fig2 = figure(2);
allvals = tracedata(traceid,:,6);
allvals = smoothignorenans(allvals,10);
time = 1:414;
plot(time,allvals)
hline = line('XData',time(1), 'YData',allvals(1), 'Color','r', ...
    'Marker','o', 'MarkerSize',6, 'LineWidth',2);

for i = 1:414
    set(hline,'XData',time(i),'Ydata',allvals(i));
    drawnow 
    pause(delay)
    f(i) = getframe;
end

