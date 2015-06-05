function output_txt = datatip_image_singletrace_angie(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

global path moviename allx ally allxorg allyorg plotfignum
nucr = 24;           %average radius of cell's nucleus
cropmult = 10;      %how far around the cell to crop the image
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
rawimage = single(imread([path,'CroppedProcessed\',moviename,'nucedge_',num2str(frame),'.tif']));
rawimage2 = single(imread([path,'CroppedProcessed\',moviename,'_pcna_',num2str(frame),'.tif']));
%rawimage = single(imread([path,'CroppedProcessed\',moviename,'\nucedge_',num2str(frame),'.tif']));
%rawimage2 = single(imread([path,'CroppedProcessed\',moviename,'\nuc_',num2str(frame),'.tif']));

%%%%% crop image %%%%%%%%%%
[m,n]=size(rawimage);
cx = int16(ally(str2num(traceid),frame));   %bestsp stores x,y data oppositely
cy = int16(allx(str2num(traceid),frame));
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
cellsinimage = find(allxorg(:,frame)>miny & allxorg(:,frame)<maxy & allyorg(:,frame)>minx & allyorg(:,frame)<maxx);
for i=cellsinimage
    text(double(int16(allxorg(i,frame))-miny)*imagescale,double(int16(allyorg(i,frame))-minx)*imagescale,num2str(i),'horizontalalignment','center','color','r','fontsize',10,'fontweight','bold');
end
%%%%% return control to plot %%%%%
figure(plotfignum);            %returns control to plot