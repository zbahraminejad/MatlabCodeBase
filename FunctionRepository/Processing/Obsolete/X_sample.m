cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
SF=1;EF=219; %Sabrina20x:208 Steve20x:219
nucr=16; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
movie=1;
%%% setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,shot,'_alldata'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%% sample specification %%%%%%%%%%%
row=2;
col=1;
site=1;
track=[779];
%%%%%  attributes  %%%%%%%%
mon=zeros(120,1);
msn=zeros(120,1);
sbn=zeros(120,1);
ssn=zeros(120,1);
aor=zeros(120,1);
anr=zeros(120,1);
sor=zeros(120,1);
snr=zeros(120,1);
%sc=zeros(120,1);
%mctd=zeros(120,1);
%sctd=zeros(120,1);
%areacell=zeros(120,1);
areabignuc=zeros(120,1);
areasmallnuc=zeros(120,1);
%%%%%%%%%%%%%%%%%%%%%%

%% setup tempwell
shottrack=[shot,'_track',num2str(track)];
if movie
    M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
    M.FrameRate = 4;
    open(M);
end
tempread=imread([cpdir,shot,'_nucedge_',num2str(1),'.tif']);
[m,n]=size(tempread);
for f=best_rc(track,1):best_rc(track,3)
    %% read in pre-cropped images
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    %CEs_or=single(imread([tifdir,'/',movieName,'_geminin_',num2str(f),'.tif']));
    %%%%%  focus on track   %%%%%%%%%%%%%%
    cx=int16(bestsp{f}(track,2));
    cy=int16(bestsp{f}(track,1));
    miny=cy-7*nucr; maxy=cy+8*nucr;
    minx=cx-7*nucr; maxx=cx+8*nucr;
    if minx<1
        minx=1; maxx=1+16*nucr;
    end
    if miny<1
        miny=1; maxy=1+16*nucr;
    end
    if maxx>m
        maxx=m; minx=m-16*nucr;
    end
    if maxy>n
        maxy=n; miny=n-16*nucr;
    end
    DAs_or=DAs_or(minx:maxx,miny:maxy);
    REs_or=REs_or(minx:maxx,miny:maxy);
    %CEs_or=CEs_or(minx:maxx,miny:maxy);
    [height,width]=size(DAs_or);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% image processing
    DAs_bl=log(imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric')); 
    DAs_bs=bgsub(DAs_bl,10*nucr,0.05);  
    REs_bl=(imfilter(REs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    REs_bs=bgsub(REs_bl,10*nucr,0.05);
    %CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    %CEs_bs=bgsub(CEs_bl,10*nucr,0.05);
    
    %% get data
    DAs_ma=getdapimask(DAs_bs,nucr);  %make the mask. this segmentation function is the engine of Feng Chiao's script
    DAs_la=bwlabel(DAs_ma);  %labels the objects with numbers

    centralnuc=DAs_la==DAs_la(floor(height/2),floor(width/2));
    DAs_la=bwlabel(centralnuc);
    
    orgring=imdilate(DAs_la,strel('disk',4,8)) - imdilate(DAs_la,strel('disk',1,8));
    BAs_la = imdilate(DAs_la,strel('disk',2,8));
    MAs_la=imerode(DAs_la,strel('disk',4,8));
    %newring=DAs_la-imerode(DAs_la,strel('disk',3,8)); %corrected to match width of orgring
    newring=BAs_la-imerode(BAs_la,strel('disk',3,8)); %corrected to match width of orgring
    
    orgmoviering=imdilate(MAs_la,strel('disk',1,8))-MAs_la;
    orgmoviering=(orgmoviering>0)*10;
    modmoviering=imdilate(MAs_la,strel('disk',1,8))-MAs_la;
    modmoviering=(modmoviering>0)*10;
    
    
    DAs_da=regionprops(DAs_la,'PixelIdxList');
    BAs_da=regionprops(BAs_la,'PixelIdxList','Area');
    MAs_da=regionprops(MAs_la,'PixelIdxList','Area');
    orgringxypos=regionprops(orgring, 'PixelIdxList');
    newringxypos=regionprops(newring, 'PixelIdxList'); 
    
    %{
    if ismember(f,[17:20])
        REs_ma=getangiemask_tf_intensityseg(REs_bs,DAs_la);
    else
        REs_ma=getangiemask_tf_intensity(REs_bs,DAs_la);
    end
    
    REs_la=bwlabel(REs_ma);
    cytoedge=bwmorph(REs_ma,'remove');
    REs_da=regionprops(REs_la,'PixelIdxList','Area');
    %}
    
    %% Summary Stats
    medorgnuc=median(REs_or(DAs_da(1).PixelIdxList));
    
    %medsmallnuc=median(REs_or(MAs_da(1).PixelIdxList));
    medsmallnuc = prctile(REs_or(MAs_da(1).PixelIdxList),25);
    %smallnuc = REs_or(MAs_da(1).PixelIdxList);
    %bottomsmallnuc = smallnuc(smallnuc<=prctile(smallnuc,50));
    %medsmallnuc = mean(bottomsmallnuc);
    
    sumbignuc=sum(sum(REs_or(BAs_da(1).PixelIdxList)));
    sumsmallnuc=sum(sum(REs_or(MAs_da(1).PixelIdxList)));
    
    allringpixels=REs_or(orgringxypos(1).PixelIdxList);
    topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
    avgorgring=mean(topringpixels);  %take mean of pixels in the top 50th pctile
    sumorgring=sum(sum(allringpixels));
    
    allringpixels=REs_or(newringxypos(1).PixelIdxList);
    topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
    avgnewring=mean(topringpixels);  %take mean of pixels in the top 50th pctile
    sumnewring=sum(sum(allringpixels));
    
    %sumcyto=sum(sum(REs_or(REs_da(1).PixelIdxList)));
    %medctd=median(CEs_or(DAs_da(1).PixelIdxList));
    %sumctd=sum(sum(CEs_or(DAs_da(1).PixelIdxList)));
    
    mon(f)=medorgnuc;
    msn(f)=medsmallnuc;
    sbn(f)=sumbignuc;
    ssn(f)=sumsmallnuc;
    aor(f)=avgorgring;
    anr(f)=avgnewring;
    sor(f)=sumorgring;
    snr(f)=sumnewring;
    %sc(f)=sumcyto;
    %mctd(f)=medctd;
    %sctd(f)=sumctd;
    %areacell(f)=REs_da(1).Area;
    areabignuc(f)=BAs_da(1).Area;
    areasmallnuc(f)=MAs_da(1).Area;
    
    
    %% Add frame to movie
    if movie
        %tempframe=imadjust(mat2gray(orgmoviering));  %red; to make a movie incluing nuclear marker
        tempframe=imadjust(mat2gray(DAs_or));
        
        tempframe(:,:,2)=imadjust(mat2gray(REs_or));
        %tempframe(:,:,2)=0;
        
        %tempframe(:,:,3)=imadjust(mat2gray(cytoedge));
        tempframe(:,:,3)=0;
        
        writeVideo(M,im2frame(tempframe));
    end
end

%%
save([datadir,shottrack,'attributes.mat'],'mon','msn','sbn','ssn','aor','anr','sor','snr','areabignuc','areasmallnuc');
if movie
    close(M);
end
