function output_txt = datatip_image_channelx(~,event_obj)

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel channelnames framesperhr
winrad=5*nucr;    %window radius: default 5
imagescale=4;     %resize image so it's bigger
% framesperhr=5;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagefignum=plotfignum+1;
pos=get(event_obj,'Position');
frame=pos(1)*framesperhr;
if frame==immunoframe
    framestring='stain';
else
    framestring=num2str(frame);
end
signal=num2str(pos(2),4);
output_txt={['X: ',framestring],['Y: ',signal]};
%%% get x and y coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tar=get(event_obj,'Target');
traceid=str2double(get(tar,'DisplayName'));  %returns cellid
allx=tracedata(:,frame,1)-jitters(frame,1); ally=tracedata(:,frame,2)-jitters(frame,2);
cx=int16(allx(traceid));
cy=int16(ally(traceid));
%%%%% load image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(imagefignum);  %makes image figure the current figure
clf;
switch channel
    case 1
        %filecode='CFP_';
        filecode=channelnames{1};
        maxval=max(tracedata(traceid,:,5),[],2);
        bgperctile=50;
    case 2
        %filecode='YFP_';
        filecode=channelnames{2};
        maxval=max(tracedata(traceid,:,6),[],2);
        bgperctile=10;
    case 3
        %filecode='RFP_';
        filecode=channelnames{3};
        maxval=max(tracedata(traceid,:,7),[],2);
        bgperctile=50;
end
raw=single(imread([rawdir,filecode,framestring,'.tif']));
%mask=single(imread([maskdir,'nucedge_',framestring,'.tif']));
%%%%% crop image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw);
% minpix=median(prctile(raw,0));  %sets absolute intensities. default 1
% maxpix=median(prctile(raw,25)); %default 99
minx=max([cx-winrad 1]);
maxx=min([cx+winrad width]);
miny=max([cy-winrad 1]);
maxy=min([cy+winrad height]);
%cropmask=mask(miny:maxy,minx:maxx);
cropraw=raw(miny:maxy,minx:maxx);
%%% calc background %%%
% bgmask=imfill(cropmask,'holes');
% bgraw=cropraw.*~bgmask;
% bgraw(bgraw==0)=NaN;
% cropbg=prctile(bgraw(:),bgperctile);
% cropraw=cropraw-prctile(cropraw(:),10);
% cropraw=cropraw-cropbg;
%cropraw(cropraw>maxval)=maxval;
%cropmask=mat2gray(imresize(cropmask,imagescale));
% cropraw=mat2gray(imresize(cropraw,imagescale),[double(minpix) double(maxpix)]);
% cropraw=imresize(mat2gray(cropraw),imagescale);
cropmask=mat2gray(imresize(cropraw,imagescale),[0 150]);
cropmask(:,:,2)=0;
cropmask(:,:,3)=0;
image(cropmask);
%%%%% add cellid over cells in cropped image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsinimage=find(ally>miny & ally<maxy & allx>minx & allx<maxx);
for i=cellsinimage
    %text(double(int16(allx(i))-minx)*imagescale,double(int16(ally(i))-miny)*imagescale,num2str(gatedtraces(i,frame,end-1)),'horizontalalignment','center','color','r','fontsize',10,'fontweight','bold');
end
%%%%% return control to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis off
figure(plotfignum);            %returns control to plot