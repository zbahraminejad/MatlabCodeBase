function A_crop(row,col)
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_10x_Steve\';  %folder containing the movies folders [CHANGE]
%path = 'h:\Documents\Timelapse\Timescape\20120807_12drugs\';
path = 'h:\Documents\Timelapse\Timescape\Steve_siCyclin\';
%path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%path = 'h:\Documents\Timelapse\Timescape\20130207_Sabrina-MCF10A-p21K\';
%datadir = ([path,'Data\']);
rawdir = ([path,'Raw\']);
cpdir = [path,'CroppedProcessed\'];
SF=1;EF=140; %Sabrina20x:208 Steve20x:219 Steve10x:111 Steve&Sabrina: 240
DAroot = '_CFP_';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
REroot = '_YFP_';  %hDHB Pilot:'EYFP' Steve20x:'EYFP' Steve10x:'EYFP' Steve&Sabrina:'YFP'
CEroot = '_mRFP_';%gemenin/cdt1 Pilot:'ECFP' Steve20x:'mRFP1' Steve10x:'mRFP1' Steve&Sabrina:'CFP'
DAname = '_nuc_';
%DAedgename = '_nucedge_';
REname = '_hDHB_';
CEname = '_geminin_';  %Pilot:'cdt1' Steve20x:'geminin' Steve10x:'geminin' Steve&Sabrina:'cdt1'

%%% setup tempwell
shot=[num2str(row),'_', num2str(col)];

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
    
    %%% patch for blank image %%%%%%%%%%%%%%
    %{
    g=f;
    if f==82    %Steve20x:82 Steve&Sabrina:41
        continue
    elseif f>82
        g=f-1;
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% reading images
    DArawfile=[shot,DAroot,num2str(f),'.tif'];
    RErawfile=[shot,REroot,num2str(f),'.tif'];
    CErawfile=[shot,CEroot,num2str(f),'.tif'];    
    DAwritefile=[shot,DAname,num2str(f),'.tif'];  %f-->g if using patch
    REwritefile=[shot,REname,num2str(f),'.tif'];  %f-->g if using patch
    CEwritefile=[shot,CEname,num2str(f),'.tif'];  %f-->g if using patch
    
    DAs_or=single(imread([rawdir,DArawfile]));
    REs_or=single(imread([rawdir,RErawfile]));
    CEs_or=single(imread([rawdir,CErawfile])); 

    %%% cropping.  
    DAs_or=CropJitter(DAs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));
    REs_or=CropJitter(REs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));
    CEs_or=CropJitter(CEs_or, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));

    %%% save cropped and processed images             
    imwrite(uint16(DAs_or),[cpdir,DAwritefile]);
    imwrite(uint16(REs_or),[cpdir,REwritefile]);
    imwrite(uint16(CEs_or),[cpdir,CEwritefile]);            

    toc(time2)
end
cd ..\Processing; %return to this directory