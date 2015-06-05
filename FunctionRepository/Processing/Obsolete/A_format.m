cd ..\Functions; %change directory for function calls
path = 'h:\Documents\CDK4 Sensor\20130715\';
%camdir = [path,'Camera\2013-04-26\1300\'];
camdir = 'E:\External JY\MCF10A\mChy-p21\2013-07-14\1682\';
rawdir = [path,'Raw\'];
w1root = 'CFP';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
w2root = 'YFP';  %hDHB Pilot:'EYFP' Steve20x:'EYFP' Steve10x:'EYFP' Steve&Sabrina:'YFP'
w3root = 'TexasRed';%gemenin/cdt1 Pilot:'ECFP' Steve20x:'mRFP1' Steve10x:'mRFP1' Steve&Sabrina:'CFP'

times = 1:117;
%rows = {'B' 'C' 'D' 'E' 'F' 'G' 'H'};
rows = {'C'};
cols = {'05'};
site = '1';

for rowidx=1:length(rows)
    row = rows{rowidx};
    for colidx=1:length(cols)
        col = cols{colidx};
        shot = [row,'_',col,'_',site];
        system(['mkdir ',rawdir,shot]);
        shotdir = [rawdir,shot,'\'];
    	for time=times
            timedir = [camdir,'TimePoint_',num2str(time),'\'];
            %system(['del ',timepath,'*thumb* /Q']);
            w1old = ['*_',row,col,'_w1*'];
            w2old = ['*_',row,col,'_w2*'];
            w3old = ['*_',row,col,'_w3*'];
            w1new = [w1root,'_',num2str(time),'.tif'];
            w2new = [w2root,'_',num2str(time),'.tif'];
            w3new = [w3root,'_',num2str(time),'.tif'];
            %system(['copy ',timedir,w1old,' ',shotdir,w1new]);
            %system(['copy ',timedir,w2old,' ',shotdir,w2new]);
            %system(['copy ',timedir,w3old,' ',shotdir,w3new]);
            system(['move ',timedir,w1old,' ',rawdir,shot]);
            system(['move ',timedir,w2old,' ',rawdir,shot]);
            system(['move ',timedir,w3old,' ',rawdir,shot]);
            system(['ren ',shotdir,w1old,' ',w1new]);
            system(['ren ',shotdir,w2old,' ',w2new]);
            system(['ren ',shotdir,w3old,' ',w3new]);
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