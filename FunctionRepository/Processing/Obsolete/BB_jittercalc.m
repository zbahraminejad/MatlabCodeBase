%function BB_jittercalc(row,col,site)
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls

%row=3;col=3;site=1;
row='B';col='03';site='1';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
SF=1;EF=230;
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
datadir = [path,'Data_Temp\'];
directory = 1;
if directory==1
    cpdir = [path,'Processed\',shot,'\'];
else
    cpdir = [path,'Processed'];
end

mask = 'nucedge&ring_';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
x=zeros(1,EF); y=zeros(1,EF);
for f=SF:EF
    im1=single(imread([cpdir,mask, num2str(f), '.tif']));
    im1=imfill(im1);
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
        %fprintf([shot,' frame %0.0f\n'],f);
        %fprintf('offset = %0.0f,%0.0f\n',xx,yy);
    end
    im0=im1;
end
save([datadir,'jitter_', shot, '.mat'],'x','y');

cd([codepath,'Processing\']); %return to this directory