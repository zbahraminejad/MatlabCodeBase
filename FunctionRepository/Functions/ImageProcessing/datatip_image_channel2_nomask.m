function output_txt = datatip_image_channel2(~,event_obj)

global rawdir maskdir gatedtraces jitters plotfignum immunoframe nucr
winrad=5*nucr;    %window radius
imagescale=4;     %resize image so it's bigger
framesperhr=5;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagefignum=plotfignum+1;
pos=get(event_obj,'Position');
frame=pos(1)*framesperhr;
if frame==immunoframe
    framestring='poststain';
else
    framestring=num2str(frame);
end
signal=num2str(pos(2),4);
output_txt={['X: ',framestring],['Y: ',signal]};
%%% get x and y coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tar=get(event_obj,'Target');
traceid=str2double(get(tar,'DisplayName'));  %returns cellid
allx=gatedtraces(:,frame,1)-jitters(frame,1); ally=gatedtraces(:,frame,2)-jitters(frame,2);
cx=int16(allx(traceid));
cy=int16(ally(traceid));
%%%%% load image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(imagefignum);  %makes image figure the current figure
clf;
%mask=single(imread([maskdir,'nucedge_',framestring,'.tif']));
raw=single(imread([rawdir,'YFP_',framestring,'.tif']));
%raw=single(imread([rawdir,'ECFP_',framestring,'.tif']));
%%%%% crop image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw);
minpix=median(prctile(raw,1));  %sets absolute intensities
maxpix=median(prctile(raw,99));
minx=max([cx-winrad 1]);
maxx=min([cx+winrad width]);
miny=max([cy-winrad 1]);
maxy=min([cy+winrad height]);
%cropmask=mask(miny:maxy,minx:maxx);
cropraw=raw(miny:maxy,minx:maxx);
%cropmask=imadjust(mat2gray(imresize(cropmask,imagescale)));
cropraw=mat2gray(imresize(cropraw,imagescale),[double(minpix) double(maxpix)]);
cropmask=cropraw;
cropmask(:,:,2)=0;
%cropmask(:,:,2)=cropraw;
cropmask(:,:,3)=0;
image(cropmask);
%%%%% add cellid over cells in cropped image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellsinimage=find(ally>miny & ally<maxy & allx>minx & allx<maxx);
for i=cellsinimage
    text(double(int16(allx(i))-minx)*imagescale,double(int16(ally(i))-miny)*imagescale,num2str(gatedtraces(i,frame,end-1)),'horizontalalignment','center','color','r','fontsize',10,'fontweight','bold');
end
%%%%% return control to plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(plotfignum);            %returns control to plot