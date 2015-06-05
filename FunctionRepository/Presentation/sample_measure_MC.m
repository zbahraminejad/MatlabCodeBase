function sample_measure_MC()
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
%%% set average nuclear radius %%%%%%%%%%%%%%%%%%%%%%%
%OldAxon: MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%IX-Micro: MCF10A/10x:12 20x:25
nucr=12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movie=1;
%%%% Sample specification %%%%%%%%%%%
row=1;
col=11;
site=1;
track=19;
%%% Setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
shottrack=[shot,'_track',num2str(track)];
load([datadir,shot,'_alldata_HS_OP'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
colorcode = jet(128);
colorcode1 = colorcode(:,1);
colorcode2 = colorcode(:,2);
colorcode3 = colorcode(:,3);
crop = 0;   %0: no cropping, 1: crop
windowsize = 10;
%windowsize = 30;
%windowsize = 10;
windowhalf = windowsize/2;
cropsize = 1;
crophalf = cropsize/2;
difhalf = windowhalf-crophalf;
framezero=63;   %Steve10xGeminin:63

%%% Set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[height,width]=size(tempread);
blocknum=3;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
minnucarea=pi*(nucr/3)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if movie
    M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
    M.FrameRate = 4;
    open(M);
end

%for f=best_rc(track,1)+108:best_rc(track,3)
for f=220
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
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_bl=log(DAs_or);
    %DAs_filt=imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric');
    %DAs_bsl=bgsub_MC(DAs_bl,blockheight,blockwidth);
    %DAs_bsf=bgsub_MC(DAs_filt,blockheight,blockwidth);
    DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    %DAs_bs=bgsub_MC(DAs_or,blockheight,blockwidth);
    
    %{
    cytogradients=edge(DAs_bs,'canny',[0.0063 0.0256]); %[.0063 .0156]
    cytogradients=bwmorph(cytogradients,'diag'); %removes 8-connectivity of background
    cytogradients=imfill(cytogradients,'holes');
    cytogradients=bwmorph(cytogradients,'thin',Inf);
    %}
    
    REs_bs=bgsub_MC(log(REs_or),blockheight,blockwidth);
    %%% Extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DAs_pad=getdapimask(DAs_bs,nucr);          %FCT algorithm
    %DAs_pad=getnucmask(DAs_bs,nucr);           %MC concave detector
    DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector 
    %DAs_pad=product_getnucmask_rev03(DAs_bs,nucr); 
    
    [DAs_da,realnuc_la,finalcytoring]=buildcytoring(DAs_pad,REs_bs,nucr);
    fcr_da=regionprops(finalcytoring,'PixelIdxList');
    
    %%% Screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    numcells=size(DAs_da,1);
    nucexclude=zeros(numcells,1);
    for k=1:numcells
        if DAs_da(k).Area < minnucarea
            nucexclude(k)=1;
        end
    end
    nucexclude=find(nucexclude);
    DAs_da(nucexclude)=[];
    fcr_da(nucexclude)=[];
    %}
    %%% Only for visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    realnuc_la(ismember(realnuc_la,nucexclude))=0;
    finalcytoring(ismember(finalcytoring,nucexclude))=0;
    nucmoviering=bwmorph(realnuc_la,'remove');
    cytomoviering=logical(finalcytoring);
    %%% Compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    mz = size(DAs_da,1);
    XX=zeros(mz,1);YY=zeros(mz,1);  %define vectors which will hold the info from the structured array, DAs_da 
    AC=zeros(mz,1);%PP=zeros(size(DAs_da,1),1);DI=zeros(size(DAs_da,1),1);
    DD=zeros(mz,1);RR=zeros(mz,1); %FR=zeros(size(DAs_da,1),1); %DD is a vector of median dapi intensities for each object in log scale
    avgringyfp=zeros(mz,1);
    for cc=1:mz   %run a loop to put the info from the structured array in to the vectors; this is for each cell
        XX(cc,1)=DAs_da(cc).Centroid(1);  %x value of centroid
        YY(cc,1)=DAs_da(cc).Centroid(2);  %y value of centroid
        AC(cc,1)=DAs_da(cc).Area;
        DD(cc,1)=median(DAs_bs(DAs_da(cc).PixelIdxList));
        RR(cc,1)=median(REs_bs(DAs_da(cc).PixelIdxList));  %nuc yfp intensity

        allringpixels=REs_bs(ringxypos(cc).PixelIdxList);
        topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
        avgringyfp(cc,1)=mean(topringpixels);  %take mean of pixels in the top 50th pctile
    end
    %}
    
    %%% Add frame to movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if movie
        
        tempframe=imadjust(mat2gray(DAs_bs));
        tempframe(:,:,2)=nucmoviering;
        %tempframe(:,:,3)=imadjust(mat2gray(nucmoviering));
        %tempframe(:,:,3)=nucmoviering;
        tempframe(:,:,3)=0;
        %}
        
        
        %{
        heatimage = imadjust(mat2gray(REs_bs));
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
        
        %tempframe = tempframe(1+difhalf*nucr:end-difhalf*nucr,1+difhalf*nucr:end-difhalf*nucr,:);
        %imshow(imresize(tempframe,4));
        %imshow(tempframe);
        writeVideo(M,im2frame(tempframe));
    end
    toc(time1)
end

%save([datadir,shottrack,'attributes.mat'],'mon','msn','sbn','ssn','aor','anr','sor','snr','areabignuc','areasmallnuc');
if movie
    close(M);
    %set(gcf,'PaperPosition',[0 0 20 20]);
    %saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory

end