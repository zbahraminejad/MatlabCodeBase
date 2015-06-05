row=4;col=6;site=14;track=26;
%row=5;col=7;site=11;track=99;
%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%% set filepaths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagepath='E:\';
imagepath='D:\Images\';
projectpath='D:\Documents\Projects\';
experimentpath='20130607 p21dCy2\20140617 PCNAvsp21dCy1dK\';

maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
datadir=([projectpath,experimentpath,'Data\']);
namemask='nucedge_';
name1='H2B_';
name2='DHB_';
name3='p21dCy1dK_';
%name1='CFP_real_';
%name2='YFP_real_';
%name3='TexasRed_real_';
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy);
%%% set viewing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=25; %default 50
imagescale=5;
nucr=4;
compression=4;
framesperhr=5;
channel2=6;
channel3=7;
channel2log=6;
linewidthvar=3;
fontsizevar=20;

shottrack=[shot,'_track',num2str(track)];
M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
M.FrameRate = 4;
open(M);
tempread=imread([maskdir,namemask,num2str(1),'.tif']);
[height,width]=size(tempread);
numframes=tracestats(track,2);
nanvec=ones(numframes,1)*NaN;
firstframe=tracestats(track,1);
lastframe=tracestats(track,2);
t2=tracedata(track,firstframe:lastframe,channel2);
min2=200; max2=prctile(t2,90);
t2log=tracedata(track,firstframe:lastframe,channel2log);
min2log=min(t2log); max2log=max(t2log);
t3=tracedata(track,firstframe:lastframe,channel3);
%min3=200; max3=max(t3);
min3=0; max3=max(t3);
finalframesize=(1+2*winrad)*imagescale;
for f=firstframe:lastframe
    %%% read in pre-cropped images %%%%%%%
    nucmask=single(imread([maskdir,namemask,num2str(f),'.tif']));
    %raw1=single(imread([rawdir,name1,num2str(f),'.tif']));
    %real2=single(imread([rawdir,name2,num2str(f),'.tif']));
    %real3=single(imread([rawdir,name3,num2str(f),'.tif']));
    real2=single(imread([maskdir,name2,num2str(f),'.tif']));
    real3=single(imread([maskdir,name3,num2str(f),'.tif']));
    %%%%%  focus on track   %%%%%%%%%%%%%%
    cx=int16(tracedata(track,f,1)-jitters(f,1));
    cy=int16(tracedata(track,f,2)-jitters(f,2));
    miny=cy-winrad; maxy=cy+winrad;
    minx=cx-winrad; maxx=cx+winrad;
    if minx<1
        minx=1; maxx=1+2*winrad;
    end
    if miny<1
        miny=1; maxy=1+2*winrad;
    end
    if maxx>width
        maxx=width; minx=width-2*winrad;
    end
    if maxy>height
        maxy=height; miny=height-2*winrad;
    end
    nucmask=nucmask(miny:maxy,minx:maxx);
    %raw1=raw1(miny:maxy,minx:maxx);
    index=f-firstframe+1;
    real2=real2(miny:maxy,minx:maxx);
    real3=real3(miny:maxy,minx:maxx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Add frame to movie
    %tempframe=mat2gray(imresize([nucmaskorg],imagescale));
    %tempframe(:,:,2)=mat2gray(imresize(real2,imagescale),range2);
    %tempframe(:,:,3)=0;
    %tempframe(:,:,1)=mat2gray(imresize(nucmask,imagescale));
    
    x=firstframe:f;
    xtime=x/framesperhr;
    y2=tracedata(track,x,channel2log);
    y3=tracedata(track,x,channel3);
    
    emptyframe=zeros(finalframesize,finalframesize);
    frame2=mat2gray(imresize(real2,imagescale),[min2 max2*0.9]);
    %plot(xtime,y2); xlim([0 lastframe/framesperhr]); ylim([min2 max2]);
    plot(xtime,y2,'linewidth',linewidthvar); xlim([0 lastframe/framesperhr]); ylim([min2log max2log]);
    xlabel('Time (hr)'); ylabel('RFU');
    hxl=get(gca,'xlabel'); set(hxl,'fontsize',fontsizevar);
    hyl=get(gca,'ylabel'); set(hyl,'fontsize',fontsizevar);
    set(gcf,'color','w');
    set(gca,'fontsize',fontsizevar);
    
    plot2=getframe(gcf); close(gcf);
    plot2=plot2.cdata; plot2=squeeze(plot2(:,:,1));
    plot2=mat2gray(imresize(plot2,[finalframesize finalframesize]));
    
    frame3=mat2gray(imresize(real3,imagescale),[min3 max3*0.9]);
    plot(xtime,y3,'linewidth',linewidthvar); xlim([0 lastframe/framesperhr]); ylim([min3 max3]);
    xlabel('Time (hr)'); ylabel('RFU');
    hxl=get(gca,'xlabel'); set(hxl,'fontsize',fontsizevar);
    hyl=get(gca,'ylabel'); set(hyl,'fontsize',fontsizevar);
    set(gcf,'color','w');
    set(gca,'fontsize',fontsizevar);
    plot3=getframe(gcf); close(gcf);
    plot3=plot3.cdata; plot3=squeeze(plot3(:,:,1));
    plot3=mat2gray(imresize(plot3,[finalframesize finalframesize]));
    
    framemask=mat2gray(imresize(nucmask,imagescale));
    
%     tempframe=mat2gray(imresize([nucmask,nucmask;nucmask,nucmask],imagescale));
%     tempframe(:,:,2)=[frame2,plot2;frame3,plot3];
%     tempframe(:,:,3)=[emptyframe,plot2;emptyframe,plot3];
%     tempframe(:,:,1)=[framemask,plot2;framemask,plot3];

    tempframe=mat2gray(imresize([nucmask,nucmask;nucmask,nucmask],imagescale));
    tempframe(:,:,2)=[frame2,plot2;emptyframe,plot3];
    tempframe(:,:,3)=[framemask,plot2;framemask,plot3];
    tempframe(:,:,1)=[emptyframe,plot2;frame3,plot3];
    
    writeVideo(M,im2frame(tempframe));
end
close(M);
fprintf('min2 = %0.0f\n',min(min2));
fprintf('max2 = %0.0f\n',max(max2));
fprintf('min3 = %0.0f\n',min(min3));
fprintf('max3 = %0.0f\n',max(max3));
