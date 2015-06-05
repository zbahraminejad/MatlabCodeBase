function makepanel(coor,rawimg,mask,maxnum,tileradius,satthresh)
[height,width]=size(rawimg);
numcells=size(coor,1);
if numcells>=maxnum
    numcells=maxnum;
end
sidenum=ceil(sqrt(numcells));
unitlength=2*tileradius+1;
sidelength=sidenum*unitlength;
panel_raw=zeros(sidelength,sidelength);
panel_nucedge=zeros(sidelength,sidelength);
panel_margin=zeros(sidelength,sidelength);
panel_margin([1 end],:)=1; panel_margin(:,[1 end])=1;
for i=1:numcells
    %xcoor=floor(nuc_rp(cellset(i)).Centroid(1));
    %ycoor=floor(nuc_rp(cellset(i)).Centroid(2));
    xcoor=floor(coor(i,1));
    ycoor=floor(coor(i,2));
    miny=ycoor-tileradius; maxy=ycoor+tileradius;
    minx=xcoor-tileradius; maxx=xcoor+tileradius;
    if minx<1
        minx=1; maxx=1+2*tileradius;
    end
    if miny<1
        miny=1; maxy=1+2*tileradius;
    end
    if maxx>width
        maxx=width; minx=width-2*tileradius;
    end
    if maxy>height
        maxy=height; miny=height-2*tileradius;
    end
    tile_raw=rawimg(miny:maxy,minx:maxx);
    tile_nucedge=mask(miny:maxy,minx:maxx);
    %tile_raw=rawimg(ycoor-tileradius+1:ycoor+tileradius,xcoor-tileradius+1:xcoor+tileradius);
    %tile_nucedge=mask(ycoor-tileradius+1:ycoor+tileradius,xcoor-tileradius+1:xcoor+tileradius);
    paneltiley=ceil(i/sidenum);
    paneltilex=mod(i,sidenum);
    if paneltilex==0
        paneltilex=sidenum;
    end
    panelposx=unitlength*(paneltilex-1)+1;
    panelposy=unitlength*(paneltiley-1)+1;
    panelxmax=panelposx+unitlength-1;
    panelymax=panelposy+unitlength-1;
    panel_raw(panelposy:panelymax,panelposx:panelxmax)=tile_raw;
    panel_nucedge(panelposy:panelymax,panelposx:panelxmax)=tile_nucedge;
    panel_margin(panelposy:panelymax,[panelposx panelxmax])=1;
    panel_margin([panelposy panelymax],panelposx:panelxmax)=1;
    %panel_raw(panelposy:panelposy+unitlength-1,panelposx:panelposx+unitlength-1)=tile_raw;
    %panel_nucedge(panelposy:panelposy+unitlength-1,panelposx:panelposx+unitlength-1)=tile_nucedge;
end
%highcalibration=max(max(rawimg));  %place in upper left corner of image to calibrate for mat2gray
%lowcalibration=min(min(rawimg));   %place in lower right corner
%panel_raw(1,1)=highcalibration;
%panel_raw(end,end)=lowcalibration;
figure;
%tempframe=mat2gray(panel_raw);
sathigh=satthresh(2);
satlow=satthresh(1);
tempframe=mat2gray(panel_raw,[satlow sathigh]);
tempframe(:,:,2)=panel_nucedge;
tempframe(:,:,3)=panel_margin;
imshow(imresize(tempframe,2));
set(gca,'units','pixels');
x=get(gca,'position');
set(gcf,'units','pixels');
y=get(gcf,'position');
set(gcf,'position',[y(1) y(2) x(3) x(4)]);
set(gca,'units','normalized','position',[0 0 1 1]);
end