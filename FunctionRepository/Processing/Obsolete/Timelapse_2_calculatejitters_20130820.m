function Timelapse_2_calculatejitters(row,col,site)
%row=3;col=3;site=1;
%row='D';col='05';site='4';
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls
timejittercalc=tic;
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\';
datadir = [path,'Data\'];
cpdir = [path,'Processed\',shot,'\'];
load([datadir,'wellsss_',shot,'.mat'],'wellsss','badframes');
totalframes = size(wellsss,3);
frameidx = find(badframes==1);
EF = totalframes-length(frameidx);

mask = 'nucedge&ring_';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
x = zeros(1,EF); y = zeros(1,EF);
g = 1;
for f=1:EF
    g=g+ismember(f,frameidx);
    im1 = single(imread([cpdir,mask, num2str(g), '.tif']));
    im1 = imfill(im1);
    if mod(f,20)==0     %just to monitor progress
        fprintf([shot,' frame %0.0f\n'],f);
    end
    if f>1
        [xx,yy] = get_alignment_shift(im0,im1,100); %max shift it will search
        x(f)=x(f-1)-xx;
        y(f)=y(f-1)-yy;
    end
    im0 = im1;
    g = g+1;
end
save([datadir,'jitter_', shot, '.mat'],'x','y');  %later just add this to wellsss
toc(timejittercalc);
cd([codepath,'Processing\']); %return to this directory