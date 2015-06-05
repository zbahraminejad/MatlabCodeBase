function Background_Precalc
bgrowmat=1;bgcolmat=1:2:11;bgsitemat=1:4;
bgrefshot='1_1_1';
realshot='8_1_3';
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
projectpath='H:\Documents\Projects\';
%experimentpath='2013-12-12_BackgroundCharacterization\';
experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
rawdir = [imagepath,experimentpath,'Raw\'];
datadir = [projectpath,experimentpath,'Data\'];
filetag={'DAPI','FITC','Cy3','Cy5'};

for i=1:length(filetag)
    refraw=single(imread([rawdir,bgrefshot,'_',filetag{i},'_stain.tif']));
    highthresh=callthresh(refraw);
    refraw=removesmears(refraw,highthresh);
    refrel=refraw(:)-nanmean(refraw(:));
    numwells=numel(bgrowmat)*numel(bgcolmat);
    meanrawrel=ones(numel(refraw),numel(bgsitemat))*NaN;
    for s=1:numel(bgsitemat)
        site=bgsitemat(s);
        rawrel=ones(numel(refraw),numwells)*NaN;
        NCS=ones(numwells,1)*NaN;        
        cc=0;
        for row=bgrowmat
            for col=bgcolmat
                cc=cc+1;
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [rawrel(:,cc),NCS(cc)]=getoffsetdiff(refrel,rawdir,shot,filetag{i},highthresh);
            end
        end
        goodwells=NCS>0.90;
        if sum(goodwells)<=2
            fprintf('Less than 3 good wells for %s\n',shot);
        end
        rawrel=rawrel(:,goodwells);
        meanrawrel(:,s)=nanmean(rawrel,2);
    end
    realraw=single(imread([rawdir,realshot,'_',filetag{i},'_stain.tif']));
    realthresh=callthresh(realraw);
    filename=[datadir,'Background_',filetag{i},'.mat'];
    save(filename,'meanrawrel','realthresh');
end

end

function [rawrel,NCS]=getoffsetdiff(refrel,rawdir,shot,filetag,highthresh)
filename=[rawdir,shot,'_',filetag,'_stain.tif'];
raw=single(imread(filename));
raw=removesmears(raw,highthresh);
rawrel=raw(:)-nanmean(raw(:));
nanidx=isnan(refrel) | isnan(rawrel);
temprefrel=refrel; temprawrel=rawrel;
temprefrel(nanidx)=[]; temprawrel(nanidx)=[];
NCS=dot(temprefrel/norm(temprefrel),temprawrel/norm(temprawrel));
end

%{
%%% debugging: view images %%%%%%%%%%
imagesc(rawimage); colorbar;
figure; imshow(mask);
%}