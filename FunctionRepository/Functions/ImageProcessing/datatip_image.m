function output_txt = datatip_image(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

global path moviename bestsp traceid plotfignum
nucr = 8;           %average radius of cell's nucleus
cropmult = 20;      %how far around the cell to crop the image
imagescale = 4;     %resize image so it's bigger

imagefignum = plotfignum+1;
%%%%% get x and y coordinates for frame & signal intensity
pos = get(event_obj,'Position');
frame = pos(1);
signal = num2str(pos(2),4);
output_txt = {['X: ',num2str(frame)],['Y: ',signal]};

%%%%% get x and y coordinates for frame & signal intensity
tar = get(event_obj,'Target');
traceid = get(tar,'DisplayName');  %returns cellid #

%%%%% load image %%%%%%%%%%
figure(imagefignum);          %makes image figure the current figure
clf;                        %clears previous image
rawimage = single(imread([path,'CroppedProcessed\',moviename,'_mRFP1_',num2str(frame),'.tif']));
rawimage2 = single(imread([path,'CroppedProcessed\',moviename,'_EYFP_',num2str(frame),'.tif']));

%%%%% crop image %%%%%%%%%%
[m,n]=size(rawimage);
bestxy = int16(bestsp{frame}(:,[1 2]));  %indices can't be decimal places
cx=bestxy(str2num(traceid),2);    %bestsp stores x,y data oppositely
cy=bestxy(str2num(traceid),1);
halfcrop = cropmult/2;
miny=cy-halfcrop*nucr; maxy=cy+halfcrop*nucr;
minx=cx-halfcrop*nucr; maxx=cx+halfcrop*nucr;
cropmult = cropmult+1;      %window is actually cropmult+1
%%%%% correct for borders %%
if minx<1
    minx=1; maxx=1+cropmult*nucr;
end
if miny<1
    miny=1; maxy=1+cropmult*nucr;
end
if maxx>m
    maxx=m; minx=m-cropmult*nucr;
end
if maxy>n
    maxy=n; miny=n-cropmult*nucr;
end
cropimage = rawimage(minx:maxx,miny:maxy);
cropimage2 = rawimage2(minx:maxx,miny:maxy);
cropimage = imadjust(mat2gray(imresize(cropimage,imagescale)));
cropimage2 = imadjust(mat2gray(imresize(cropimage2,imagescale)));
cropimage(:,:,2) = cropimage2;
cropimage(:,:,3) = 0;
image(cropimage);
%%%%% add cellid over cells in cropped image %%%%%%
cellsinimage = find(bestxy(:,1)>miny & bestxy(:,1)<maxy & bestxy(:,2)>minx & bestxy(:,2)<maxx);
for i=cellsinimage
    text(double(bestxy(i,1)-miny)*imagescale,double(bestxy(i,2)-minx)*imagescale,num2str(i),'horizontalalignment','center','color','k','fontsize',10,'fontweight','bold');
end
%%%%% return control to plot %%%%%
figure(plotfignum);            %returns control to plot