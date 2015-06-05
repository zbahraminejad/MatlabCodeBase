function reformatbydatemodified(row,col,site)
%row='B'; col='04'; site='1';
time2=tic;
cd ..\Functions; %change directory for function calls
rawdir = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\Raw\';
w1root = 'CFP';  %H2B  Pilot:'mRFP1' Steve20x:'ECFP' Steve10x:'CFP' Steve&Sabrina:'mRFP1'
w2root = 'YFP';  %hDHB Pilot:'EYFP' Steve20x:'EYFP' Steve10x:'EYFP' Steve&Sabrina:'YFP'
w3root = 'TexasRed';%gemenin/cdt1 Pilot:'ECFP' Steve20x:'mRFP1' Steve10x:'mRFP1' Steve&Sabrina:'CFP'

times = 1:149;
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
shotdir = [rawdir,shot];
listing = dir([shotdir,'\mChy*']);
totalfiles=length(times)*3;
datenumlist = zeros(totalfiles,1);
namelist = cell(totalfiles,1);
for filenum=1:totalfiles
    datenumlist(filenum) = listing(filenum).datenum;
    namelist{filenum} = listing(filenum).name;
end
[~,dateidx] = sort(datenumlist);
sortednames = namelist(dateidx);

for filenum=1:totalfiles
    time = ceil(filenum/3);
    if strcmp(sortednames{filenum}(26),'1')
        newname = [w1root,'_',num2str(time),'.tif'];
    elseif strcmp(sortednames{filenum}(26),'2')
        newname = [w2root,'_',num2str(time),'.tif'];
    elseif strcmp(sortednames{filenum}(26),'3')
        newname = [w3root,'_',num2str(time),'.tif'];
    else
        fprintf('error at file %d\n',filenum);
    end
    system(['ren ',shotdir,'\',sortednames{filenum},' ',newname]);
end
        
%%% reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system(['mkdir ',path,'Raw']);
%system(['move ',path,'testb.tif ',path,'Raw']);
%system(['ren ',path,'testf.tif testg.tif']);
%system(['copy ',path,'testa.tif ',path,'Raw\testacopied.tif']);
%system(['del ',fldir,'* /Q']);

cd ..\Processing; %return to this directory
toc(time2)