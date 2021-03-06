function A_crop_noCFP(row,col,site)
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_10x_Steve\';  %folder containing the movies folders [CHANGE]
%path = 'h:\Documents\Timelapse\Timescape\20120807_12drugs\';
%path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
path = 'h:\Documents\Timelapse\Timescape\20130207_Sabrina-MCF10A-p21KO\';
%datadir = ([path,'Data\']);
rawdir = ([path,'Raw\']);
cpdir = [path,'CroppedProcessed\'];
SF=1;EF=95; %Sabrina20x:208 Steve20x:219 Steve10x:111 Steve&Sabrina: 240
%nucr=8; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
DAroot = '_CFP_';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
REroot = '_YFP_';  %hDHB Pilot:'EYFP' Steve20x:'EYFP' Steve10x:'EYFP' Steve&Sabrina:'YFP'
DAname = '_nuc_';
REname = '_hDHB_';

%%% setup tempwell
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

%%% correcting jitters
x=zeros(1,EF); y=zeros(1,EF);

for f=SF:EF
    disp(f) %print to the screen the frame you are on
    time2=tic;

    %%% reading images
    DAs_or=single(imread([rawdir,shot, DAroot, num2str(f), '.tif']));%DA = dapi, reads image. 2D matrix of intensities for that frame

    %%% correct multiple jitters
    DI1=imadjust(mat2gray(log(DAs_or)));

    if f>SF                   
        [xx,yy]=CalcJitter(DI0,DI1);
        x(f)=x(f-1)+xx;  
        y(f)=y(f-1)+yy;
        disp([x(f), y(f)]);
    end
    DI0=DI1;
    toc(time2);
end

cropamountx=ceil(max(abs(x))+1);  %use maximum jitter to find amount to crop  
cropamounty=ceil(max(abs(y))+1);  %use maximum jitter to find amount to crop         

for f=SF:EF
    time2=tic;
    
    %%% reading images
    DArawfile=[shot,DAroot,num2str(f),'.tif'];
    RErawfile=[shot,REroot,num2str(f),'.tif'];
    DAwritefile=[shot,DAname,num2str(f),'.tif'];  %f-->g if using patch
    REwritefile=[shot,REname,num2str(f),'.tif'];  %f-->g if using patch
    
    DAs_or=single(imread([rawdir,DArawfile]));
    REs_or=single(imread([rawdir,RErawfile]));

    %%% cropping.  
    DAs_or=CropJitter(DAs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));
    REs_or=CropJitter(REs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));

    %%% save cropped and processed images             
    imwrite(uint16(DAs_or),[cpdir,DAwritefile]);
    imwrite(uint16(REs_or),[cpdir,REwritefile]);           

    toc(time2)
end
cd ..\Processing; %return to this directory