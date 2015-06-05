row = 4; col = 2; site = 1;

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
    rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'_'];
    realdir = [imagepath,experimentpath,'Real\',wellname,shot,'_'];
end

load([datadir,'tracedata_',shot,'_nolink','.mat'],'tracedata','genealogy','jitters');


%%%%%% Set subtightplot settings %%%%%%%%%%%%%%%%%%%%%%%
gapcoordinates = [0.06 0.04];
heightmargin = 0.07;
widthmargin = 0.025;

%%%%%Get x and y coordinates%%%%%
traceid=[1881]; %180 228 350 %input cell id
fig = figure(1);
set(fig,'color','white','units','normalized','outerposition',[0 0 .6 1]);
subtightplot(2,2,4,gapcoordinates,heightmargin,widthmargin)

time = 1:653;
xtime = time./5.5;
allvals1 = tracedata(traceid,:,10);
allvals1 = smoothignorenans(allvals1,6);
allvals2 = tracedata(traceid,:,11);
allvals2 = smoothignorenans(allvals2,6);
ymin1 = 0; ymax1 = 150000; ystep1 = 10000;
ymin2 = 0; ymax2 = 3; ystep2 = 1;
tracewidth = 5;
[haxes,hline1,hline2]=plotyy(xtime,allvals1,xtime,allvals2);
axes(haxes(1));
axis([xtime(1) xtime(end) ymin1 ymax1]);
set(gca,'Box','off','YAxisLocation','left','YColor','k','YTick',ymin1:ystep1:ymax1);
set(hline1,'DisplayName',num2str(i),'color',[0.9 .75 0],'linewidth',tracewidth); % was blue before 'b'

axes(haxes(2));
axis([xtime(1) xtime(end) ymin2 ymax2]);
set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2); % Ytick spacing default is --> ymin2:ystep2:ymax2
set(hline2,'DisplayName',num2str(i),'Color','r','linewidth',tracewidth);

hold on
axes(haxes(1));
hline = line('XData',time(1), 'YData',allvals(1), 'Color','r','Marker','o', 'MarkerSize',10, 'LineWidth',10);
line([113/6 113/6] ,[0 max(allvals)],'Color','k','linewidth',3,'linestyle','--');
line([385/6 385/6],[0 max(allvals)],'Color','k','linewidth',3,'linestyle','-');
hold off
set(gca,'FontName','Arial','FontSize',20);

title('Time Trace','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
xlabel('Time (hours)','FontName','Arial','FontSize',26,'Units','Normalized','Position',[.5 -0.07]);
ylabel('RFU','FontName','Arial','FontSize',26,'Units','Normalized','Position',[-.09 .5]);
text(3, 25, '+Rosi', 'Color',[0 0 0],'FontSize',30, 'HorizontalAlignment','left', 'VerticalAlignment','top');
text(48, 40, '+Insulin', 'Color',[0 0 0], 'FontSize',30, 'HorizontalAlignment','left', 'VerticalAlignment','top');

videofile ='C:\Users\teruellab\Desktop\example_pparg_cell.avi';
% M=VideoWriter(videofile,'Uncompressed AVI');
M=VideoWriter(videofile,'Motion JPEG AVI');
M.FrameRate=13;
open(M);

 %setup to get max and min intensity 
    channel = {'CFP_','YFP_'};
    realdumb=double(imread([rawdir,channel{2},num2str(200),'.tif']));
    maskdumb=double(imread([maskdir,'nucedge_',num2str(200),'.tif']));
    nucdumb = double(imread([rawdir,channel{1},num2str(200),'.tif']));
    minInt = min(realdumb(:));
    maxInt = max(realdumb(:));
for frame = 1:414;  

allx=tracedata(:,frame,1)-jitters(frame,1); ally=tracedata(:,frame,2)-jitters(frame,2);
cx=int16(allx(traceid));
cy=int16(ally(traceid));

channel = {'CFP_','YFP_'};
bgperctile=10;
maxval=max(tracedata(traceid,:,6),[],2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
real=double(imread([rawdir,channel{2},num2str(frame),'.tif']));
mask=double(imread([maskdir,'nucedge_',num2str(frame),'.tif']));
nuc = double(imread([rawdir,channel{1},num2str(frame),'.tif']));

nucr = 16;
winrad=2*nucr;    %window radius: default 5
imagescale=4;  
%%%%% crop image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(real);
%minpix=median(prctile(raw,0));  %sets absolute intensities. default 1
%maxpix=median(prctile(raw,100)); %default 99
minx=max([cx-winrad 1]);
maxx=min([cx+winrad width]);
miny=max([cy-winrad 1]);
maxy=min([cy+winrad height]);
cropmask=mask(miny:maxy,minx:maxx);
cropreal=real(miny:maxy,minx:maxx);
cropnuc = nuc(miny:maxy,minx:maxx);

% to make a bigger picture
minxLarge=max([cx-(winrad*5) 1]);
maxxLarge=min([cx+(winrad*5) width]);
minyLarge=max([cy-(winrad*5) 1]);
maxyLarge=min([cy+(winrad*5) height]);
cropmaskLarge=mask(minyLarge:maxyLarge,minxLarge:maxxLarge);
croprealLarge=real(minyLarge:maxyLarge,minxLarge:maxxLarge);

%%% calc background %%%
% bgmask=imfill(cropmask,'holes');
% bgraw=cropreal.*~bgmask;
% bgraw(bgraw==0)=NaN;
% cropbg=prctile(bgraw(:),bgperctile);
% cropraw=cropraw-prctile(cropraw(:),10);
% cropreal=cropreal-cropbg;
% cropreal(cropreal>maxval)=maxval;

% cropmask=mat2gray(imresize(cropmask,imagescale));
% cropmaskLarge=mat2gray(imresize(cropmaskLarge,imagescale));
cropmask=imresize(cropmask,imagescale);
cropmaskLarge=imresize(cropmaskLarge,imagescale);
% cropraw=mat2gray(imresize(cropraw,imagescale),[double(minpix) double(maxpix)]);
% cropraw=imresize(mat2gray(cropraw),imagescale);
cropmask(:,:,2)=mat2gray(imresize(cropreal,imagescale),[minInt maxInt]);
cropmask(:,:,3)=0;
cropmaskLarge(:,:,2)=mat2gray(imresize(croprealLarge,imagescale),[minInt maxInt]);
cropmaskLarge(:,:,3)=0;

cropnuc =  mat2gray(imresize(cropnuc,imagescale));
tempnuc = cropnuc;
cropnuc(:,:,3) = tempnuc;
cropnuc(:,:,1) = 0;
cropnuc(:,:,2) = tempnuc;
 


figure(1), subtightplot(2,2,1,gapcoordinates,heightmargin,widthmargin)
image(cropmask);
set(gca,'YTick',[])
set(gca,'XTick',[])
title('PPARg-YFP','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
% axis square

figure(1) ,subtightplot(2,2,2,gapcoordinates,heightmargin,widthmargin)
image(cropnuc)
set(gca,'YTick',[])
set(gca,'XTick',[])
title('H2B-CFP','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
% axis square

figure(1), subtightplot(2,2,3,gapcoordinates,heightmargin,widthmargin)
image(cropmaskLarge)
set(gca,'YTick',[])
set(gca,'XTick',[])
title('Group Shot','FontName','Arial','FontSize',36,'Units','Normalized','Position',[.5 1]);
% axis square

figure(1), subtightplot(2,2,4,gapcoordinates,heightmargin,widthmargin)
set(hline,'XData',time(frame)./6,'Ydata',allvals(frame));
drawnow 
axis square

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

