function FormatFiles_Immunostain(row,col,site)
row='2'; col='2'; site='1';
time2=tic;

imagepath='H:\Images\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
experimentpath='2013-11-21_IL2_pSTATpERKpAKT\Experiment_20131121\';
rawdir=[imagepath,experimentpath,'Raw\'];
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

camdir='E:\ix_micro\20131121-YT-IL2\20131121-YT-IL2-STAT5-pERKpAKT\2013-11-22\2193\';
w1root='Hoechst';
w2root='AF488';
w3root='AF546';
w4root='Cy5';
%w5root='DAPI';
stain='stain';

shotdir=[rawdir,shot];
if ~exist(shotdir,'dir')
    mkdir(shotdir);
end

filenames=getFilenames(camdir);
filenames=filenames(~boolRegExp(filenames,'thumb'));
w1old=char(filenames(boolRegExp(filenames,['_',row,col,'_s',site,'_w1'])));
w2old=char(filenames(boolRegExp(filenames,['_',row,col,'_s',site,'_w2'])));
w3old=char(filenames(boolRegExp(filenames,['_',row,col,'_s',site,'_w3'])));
w4old=char(filenames(boolRegExp(filenames,['_',row,col,'_s',site,'_w4'])));
%w5old=char(filenames(boolRegExp(filenames,['_',row,col,'_s',site,'_w5'])));    
w1new=[w1root,'_',stain,'.tif'];
w2new=[w2root,'_',stain,'.tif'];
w3new=[w3root,'_',stain,'.tif'];
w4new=[w4root,'_',stain,'.tif'];
%w5new=[w5root,'_',stain,'.tif'];    
system(['copy ',camdir,'\',w1old,' ',shotdir,'\',w1new]);
system(['copy ',camdir,'\',w2old,' ',shotdir,'\',w2new]);
system(['copy ',camdir,'\',w3old,' ',shotdir,'\',w3new]);
system(['copy ',camdir,'\',w4old,' ',shotdir,'\',w4new]);
%system(['copy ',camdir,'\',w5old,' ',shotdir,'\',w5new]);

        
%%% reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system(['mkdir ',path,'Raw']);
%system(['move ',path,'testb.tif ',path,'Raw']);
%system(['ren ',path,'testf.tif testg.tif']);
%system(['copy ',path,'testa.tif ',path,'Raw\testacopied.tif']);
%system(['del ',fldir,'* /Q']);
toc(time2)