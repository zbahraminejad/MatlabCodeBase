function Timelapse_2_calculatejitters(row,col,site)
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls
timejittercalc=tic;
%row=3;col=3;site=1;
%row='B';col='04';site='3';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
SF=1;EF=230;
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
datadir = [path,'Data\'];
cpdir = [path,'Processed\',shot,'\'];

load([datadir,'badframes.mat'],'movieindex','badframes');
wellidx=find(strcmp(movieindex,shot));
frameidx=find(badframes(wellidx,:)==1);
EF=EF-length(frameidx);

mask = 'nucedge&ring_';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
x=zeros(1,EF); y=zeros(1,EF);
g=SF;
for f=SF:EF
    if ismember(f,frameidx)
        g=g+1;
    end
    im1=single(imread([cpdir,mask, num2str(g), '.tif']));
    im1=imfill(im1);
    if mod(f,20)==0     %just to monitor progress
        fprintf([shot,' frame %0.0f\n'],f);
    end
    if f>SF
        %[xx,yy]=CalcJitter(im0,im1);
        [xx,yy]=get_alignment_shift(im0,im1,100); %max shift it will search
        if abs(xx)<200 && abs(yy)<200  %irrelevant if using get_alignment_shift
            %x(f)=x(f-1)+xx; %positive sign for CalcJitter  
            %y(f)=y(f-1)+yy;
            x(f)=x(f-1)-xx;  %reverse sign for get_alignment_shift
            y(f)=y(f-1)-yy;
        else
            fprintf(['excess jitter calc: ',shot,' frame %0.0f: %0.2f %0.2f\n'],f,xx,yy);
            x(f)=x(f-1);
            y(f)=y(f-1);
        end
        %fprintf([shot,' frame %0.0f\n'],f);
        %fprintf('offset = %0.0f,%0.0f\n',xx,yy);
    end
    im0=im1;
    g=g+1;
end
save([datadir,'jitter_', shot, '.mat'],'x','y');
toc(timejittercalc)
cd([codepath,'Processing\']); %return to this directory