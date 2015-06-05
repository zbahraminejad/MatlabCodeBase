function FormatFiles_Rename(row,col,site)
%row='B'; col='01'; site='1';
time2=tic;

imagepath='H:\Images\';
%experimentpath='2013-06-12_JY_BJ5_CycDFKBP_DHFRp21\';
%experimentpath='2013-02-22_Sabrina_MCF10Ap21null_DHFR-Chy-p21\';
experimentpath='Sabrina\';
maskdir=[imagepath,experimentpath,'Mask\'];
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%shotdir=[maskdir,shot];

w1oldroot='nucedge';
w2oldroot='EYFP';
w3oldroot='mRFP1'; %w3oldroot='"Texas Red"';
w1newroot='_nucedge';
w2newroot='YFP';
w3newroot='TexasRed';

times=1:116;
for time=times
    w1old=['_',shot,w1oldroot,'_',num2str(time),'.tif'];
    %w2old=[shot,'_',w2oldroot,'_',num2str(time),'.tif'];
    %w3old=[shot,'_',w3oldroot,'_',num2str(time),'.tif'];
    w1new=[shot,w1newroot,'_',num2str(time),'.tif'];
    %w2new=[w2newroot,'_',num2str(time),'.tif'];
    %w3new=[w3newroot,'_',num2str(time),'.tif'];
    system(['ren ',maskdir,w1old,' ',w1new]);
    %system(['ren ',maskdir,w2old,' ',w2new]);
    %system(['ren ',maskdir,w3old,' ',w3new]);
end
        
%%% reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system(['mkdir ',path,'Raw']);
%system(['move ',path,'testb.tif ',path,'Raw']);
%system(['ren ',path,'testf.tif testg.tif']);
%system(['copy ',path,'testa.tif ',path,'Raw\testacopied.tif']);
%system(['del ',fldir,'* /Q']);
toc(time2)