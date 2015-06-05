stain = 'poststain';
time2=tic;
cd ..\Functions; %change directory for function calls
rawdir = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\Raw\';
camdir = [rawdir,stain,'\'];
w1root = 'CFP';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
w2root = 'YFP';  %hDHB Pilot:'EYFP' Steve20x:'EYFP' Steve10x:'EYFP' Steve&Sabrina:'YFP'
w3root = 'TexasRed';%gemenin/cdt1 Pilot:'ECFP' Steve20x:'mRFP1' Steve10x:'mRFP1' Steve&Sabrina:'CFP'
w4root = 'Cy5';
w5root = 'DAPI';

%rows={'B'}; cols={'03'}; sites={'1'};
rows = {'B' 'C' 'D' 'E'};
cols = {'03' '04' '05' '06' '07' '08' '09' '10' '11'};
sites = {'1' '2' '3' '4'};

system(['del ',camdir,'*thumb* /Q']);
for rowidx = 1:length(rows)
    row = rows{rowidx};
    for colidx = 1:length(cols)
        col = cols{colidx};
        for siteidx = 1:length(sites)
            site = sites{siteidx};
            fprintf('%s_%s_%s\n',row,col,site);
            shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
            shotdir=[rawdir,shot];
            if ~exist(shotdir,'dir')
                mkdir(shotdir);
            end
            w1old = ['*_',row,col,'_s',site,'_w1*'];
            w2old = ['*_',row,col,'_s',site,'_w2*'];
            w3old = ['*_',row,col,'_s',site,'_w3*'];
            w4old = ['*_',row,col,'_s',site,'_w4*'];
            w5old = ['*_',row,col,'_s',site,'_w5*'];
            w1new = [w1root,'_',stain,'.tif'];
            w2new = [w2root,'_',stain,'.tif'];
            w3new = [w3root,'_',stain,'.tif'];
            w4new = [w4root,'_',stain,'.tif'];
            w5new = [w5root,'_',stain,'.tif'];
            system(['move ',camdir,w1old,' ',shotdir]);
            system(['move ',camdir,w2old,' ',shotdir]);
            system(['move ',camdir,w3old,' ',shotdir]);
            system(['move ',camdir,w4old,' ',shotdir]);
            system(['move ',camdir,w5old,' ',shotdir]);
            system(['ren ',shotdir,'\',w1old,' ',w1new]);
            system(['ren ',shotdir,'\',w2old,' ',w2new]);
            system(['ren ',shotdir,'\',w3old,' ',w3new]);
            system(['ren ',shotdir,'\',w4old,' ',w4new]);
            system(['ren ',shotdir,'\',w5old,' ',w5new]);
        end
    end
end
        
%%% reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system(['mkdir ',path,'Raw']);
%system(['move ',path,'testb.tif ',path,'Raw']);
%system(['ren ',path,'testf.tif testg.tif']);
%system(['copy ',path,'testa.tif ',path,'Raw\testacopied.tif']);
%system(['del ',fldir,'* /Q']);

cd ..\Processing; %return to this directory
toc(time2)