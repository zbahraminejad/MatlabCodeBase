function sample_measure_MC_graph()
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
%%% set average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%
%OldAxon: MCF10A/10x:8 MCF10A/20x:16
%IX-Micro: MCF10A/10x:12 MCF10A/20x:25
nucr=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movie=0;
plotdata=1;
%%%% Sample specification %%%%%%%%%%%
row=1;
col=11;
site=1;
track=214;
%track=43;
%%% Setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
shottrack=[shot,'_track',num2str(track)];
load([datadir,shot,'_alldata_HS_OP'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
colorcode = jet(128);
colorcode1 = colorcode(:,1);
colorcode2 = colorcode(:,2);
colorcode3 = colorcode(:,3);
crop = 1;   %0: no cropping, 1: crop
windowsize = 10;
%windowsize = 30;
%windowsize = 10;
windowhalf = windowsize/2;

%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[height,width]=size(tempread);
blocknum=1;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
%minnucarea=pi*(nucr/3)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if movie
    M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
    M.FrameRate = 4;
    open(M);
end

SF=best_rc(track,1);
EF=best_rc(track,3);
totalframes=EF-SF+1;
frames=SF:EF;
RR=zeros(totalframes,1);
ring=zeros(totalframes,1);

for f=SF:EF
%for f=50
    time1=tic;
    %%% Read in pre-cropped images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    %CEs_or=single(imread([tifdir,'/',movieName,'_geminin_',num2str(f),'.tif']));
    %%%%%  Focus on track   %%%%%%%%%%%%%%
    if crop ==1
        cx=int16(bestsp{f}(track,2));
        cy=int16(bestsp{f}(track,1));
        %cx=1510; cy=610;
        miny=cy-(windowhalf-1)*nucr; maxy=cy+windowhalf*nucr;
        minx=cx-(windowhalf-1)*nucr; maxx=cx+windowhalf*nucr;
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
        DAs_or=DAs_or(minx:maxx,miny:maxy);
        REs_or=REs_or(minx:maxx,miny:maxy);
        %CEs_or=CEs_or(minx:maxx,miny:maxy);
    end
    [miniheight,miniwidth]=size(DAs_or);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    REs_bs=bgsub_MC(log(REs_or),blockheight,blockwidth);
    nucbg=bgsub_MC(DAs_or,blockheight,blockwidth);
    nucbg(nucbg<1)=1;
    nucbg=log(nucbg);
    hdhbbg=bgsub_MC(REs_or,blockheight,blockwidth);
    hdhbbg(hdhbbg<1)=1;
    hdhbbg=log(hdhbbg);
    %%% Extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_pad=getdapimask(DAs_bs,nucr);          %FCT algorithm
    %DAs_pad=getnucmask(DAs_bs,nucr);           %MC concave detector
    DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector 
    %DAs_pad=product_getnucmask_rev03(DAs_bs,nucr); 
    
    [DAs_da,realnuc_la,finalcytoring]=buildcytoring_single(DAs_pad,REs_bs,nucr);
    %fcr_da=regionprops(finalcytoring,'PixelIdxList');
    
    centrallabel=realnuc_la(floor(miniheight/2),floor(miniwidth/2));
    centralnuc=realnuc_la==centrallabel;
    DAs_da=regionprops(centralnuc,'PixelIdxList');
    centralring=finalcytoring==centrallabel;
    fcr_da=regionprops(centralring,'PixelIdxList');
    
    %%% Only for visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucmoviering=bwmorph(centralnuc,'remove');
    cytomoviering=logical(centralring);

    %%% Compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RR(f)=median(REs_bs(DAs_da(1).PixelIdxList));  %nuc yfp intensity
    nucweight=DAs_bs(DAs_da(1).PixelIdxList);
    nucsig=REs_bs(DAs_da(1).PixelIdxList);
    nucsig(nucsig<0)=0;
    %sigbyweight=nucsig./nucweight;
    nucweightmax=max(nucweight);
    nucweightcoeff=nucweightmax./nucweight;
    %RR(f)=median(nucsig.*nucweightcoeff);
    RR(f)=median(REs_bs(DAs_da(1).PixelIdxList));
    ring(f)=median(REs_bs(fcr_da(1).PixelIdxList));

    %%% Add frame to movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if movie
        %{
        tempframe=imadjust(mat2gray(REs_bs));
        tempframe(:,:,2)=cytomoviering;
        %tempframe(:,:,3)=imadjust(mat2gray(nucmoviering));
        tempframe(:,:,3)=nucmoviering;
        %}
        
        nucweightcoeff=nucweightmax./DAs_bs;
        REcoeff=REs_bs.*nucweightcoeff;
        REcoeff(centralnuc==0)=0;
        
        REs_bs(centralnuc==0)=0;
        heatimage = imadjust(mat2gray(REs_bs));
        %heatimage = imadjust(mat2gray(REcoeff));
        normimage = heatimage/(max(max(heatimage)));
        codeimage = round(normimage*128);
        codeimage(codeimage==0) = 1;
        
        blackring=1;
        blackimage=nucmoviering;
        %blackimage=cytoringedges | realnuc;
        %redimage=legitedges;
        if ~blackring
            tempframe = colorcode1(codeimage);
            tempframe(:,:,2) = colorcode2(codeimage);
            tempframe(:,:,3) = colorcode3(codeimage);
        else
            tempframe = colorcode1(codeimage);
            tempframe(logical(blackimage))=0;
            %tempframe(logical(redimage))=1;
            frame2=colorcode2(codeimage);
            frame2(logical(blackimage))=0;
            %frame2(logical(redimage))=0;
            tempframe(:,:,2) = frame2;
            frame3=colorcode3(codeimage);
            frame3(logical(blackimage))=0;
            %frame3(logical(redimage))=0;
            tempframe(:,:,3) = frame3;
        end
        %}
        
        imshow(imresize(tempframe,4));
        %imshow(tempframe);
        %writeVideo(M,im2frame(tempframe));
    end
    toc(time1)
end

%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotdata
    signal1=ring./RR;
    signal2=ring;
    signal3=RR;
    plot(frames,signal1);
end

%save([datadir,shottrack,'attributes.mat'],'mon','msn','sbn','ssn','aor','anr','sor','snr','areabignuc','areasmallnuc');
if movie
    close(M);
    set(gcf,'PaperPosition',[0 0 8 4]);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory

end