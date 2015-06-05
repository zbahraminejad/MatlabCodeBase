function A_crop_directory(row,col,site)
%row='B';col='05';site='1';
SF=1;EF=230;

codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

shot=[row,'_',col,'_',site];
%path = 'h:\Documents\Timelapse\Timescape\20130424_Panel_DHB-Gem_40hr\';
%path = 'H:\Documents\Timelapse\Timescape\20130514_MCF10A-p21KO_24hr_6min\';
path = 'h:\Documents\CDK4\20130715\';
rawdir = ([path,'Raw\',shot,'\']);
cpdir = [path,'CroppedProcessed\',shot,'\'];
if ~exist(cpdir,'dir')
    mkdir(cpdir);
end

nucnameold = 'CFP_';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
DHBnameold = 'YFP_';  %hDHB Pilot:'EYFP' Steve20x:'EYFP' Steve10x:'EYFP' Steve&Sabrina:'YFP'
cdtnameold = 'TexasRed_';%gemenin/cdt1 Pilot:'ECFP' Steve20x:'mRFP1' Steve10x:'mRFP1' Steve&Sabrina:'CFP'
nucnamenew = 'nuc_';
DHBnew = 'DHB_';
cdtnamenew = 'p21_';  %Pilot:'cdt1' Steve20x:'geminin' Steve10x:'geminin' Steve&Sabrina:'cdt1'

%%% correcting jitters
x=zeros(1,EF); y=zeros(1,EF);

for f=SF:EF
    nuc_org=single(imread([rawdir,nucnameold, num2str(f), '.tif']));
    im1=imadjust(mat2gray(log(nuc_org)));
    if mod(f,20)==0     %just to monitor progress
        fprintf([shot,' frame %0.0f\n'],f);
    end
    if f>SF                   
        [xx,yy]=CalcJitter(im0,im1);
        if abs(xx)<200 && abs(yy)<200
            x(f)=x(f-1)+xx;  
            y(f)=y(f-1)+yy;
        else
            fprintf(['excess jitter calc: ',shot,' frame %0.0f: %0.2f %0.2f\n'],f,xx,yy);
            x(f)=x(f-1);
            y(f)=y(f-1);
        end
    end
    im0=im1;
end

cropamountx=ceil(max(abs(x))+1);  %use maximum jitter to find amount to crop  
cropamounty=ceil(max(abs(y))+1);  %use maximum jitter to find amount to crop         
system(['mkdir ',cpdir]);

for f=SF:EF
    fprintf('crop frame %0.0f\n',f);
    time2=tic;
    g=f;
    %%% patch for blank image %%%%%%%%%%%%%%
    %{
    if f==11 || f==41    %Steve20x:82 Steve&Sabrina:41
        continue
    elseif f>=12 && f<=40
        g=f-1;
    elseif f>=42
        g=f-2;
    end
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% reading images
    nucrawfile=[nucnameold,num2str(f),'.tif'];
    DHBrawfile=[DHBnameold,num2str(f),'.tif'];
    cdtrawfile=[cdtnameold,num2str(f),'.tif'];    
    nucwritefile=[nucnamenew,num2str(g),'.tif'];  %f-->g if using patch
    DHBwritefile=[DHBnew,num2str(g),'.tif'];  %f-->g if using patch
    cdtwritefile=[cdtnamenew,num2str(g),'.tif'];  %f-->g if using patch
    
    nuc_org=single(imread([rawdir,nucrawfile]));
    DHB_org=single(imread([rawdir,DHBrawfile]));
    cdt_org=single(imread([rawdir,cdtrawfile])); 

    %%% cropping.  
    nuc_org=CropJitter(nuc_org, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));
    DHB_org=CropJitter(DHB_org, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));
    cdt_org=CropJitter(cdt_org, cropamountx, cropamountx, cropamounty, cropamounty,  x(f), y(f));

    %%% save cropped and processed images             
    imwrite(uint16(nuc_org),[cpdir,nucwritefile]);
    imwrite(uint16(DHB_org),[cpdir,DHBwritefile]);
    imwrite(uint16(cdt_org),[cpdir,cdtwritefile]);            

    toc(time2)
end
cd([codepath,'Processing\']); %return to this directory