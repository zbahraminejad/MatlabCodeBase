function Timelapse_0_formatfiles(row,col,site)
%row='B'; col='02'; site='1';
time2=tic;

imagepath = 'H:\Images\';
experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
rawdir = [imagepath,experimentpath,'Raw\'];
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

camdir = 'E:\External_JY\MC\20130831-BJ5-p21Cy1-Geminin\2013-08-31\1904\';
w1root = 'CFP';
w2root = 'YFP';
w3root = 'TexasRed';

times = 1:12;
shotdir = [rawdir,shot];
if ~exist(shotdir,'dir')
    mkdir(shotdir);
end
for time=times
    timedir = [camdir,'TimePoint_',num2str(time),'\'];
    system(['del ',timedir,'*thumb* /Q']);
    w1old = ['*_',row,col,'_s',site,'_w1*'];
    w2old = ['*_',row,col,'_s',site,'_w2*'];
    w3old = ['*_',row,col,'_s',site,'_w3*'];
    w1new = [w1root,'_',num2str(time),'.tif'];
    w2new = [w2root,'_',num2str(time),'.tif'];
    w3new = [w3root,'_',num2str(time),'.tif'];
    system(['move ',timedir,w1old,' ',shotdir]);
    system(['move ',timedir,w2old,' ',shotdir]);
    system(['move ',timedir,w3old,' ',shotdir]);
    system(['ren ',shotdir,'\',w1old,' ',w1new]);
    system(['ren ',shotdir,'\',w2old,' ',w2new]);
    system(['ren ',shotdir,'\',w3old,' ',w3new]);
end
        
%%% reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system(['mkdir ',path,'Raw']);
%system(['move ',path,'testb.tif ',path,'Raw']);
%system(['ren ',path,'testf.tif testg.tif']);
%system(['copy ',path,'testa.tif ',path,'Raw\testacopied.tif']);
%system(['del ',fldir,'* /Q']);
toc(time2)