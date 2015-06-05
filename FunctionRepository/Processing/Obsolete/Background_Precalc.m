function Background_Precalc
%rowmat=1;colmat=[1:2:11];sitemat=1:4;
rowmat=2:4;colmat=2;sitemat=1:4;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
projectpath='H:\Documents\Projects\';
%experimentpath='2013-12-12_BackgroundCharacterization\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
experimentpath='2013-12-13_pRbAbCharacterization\Experiment_20131217\';
rawdir = [imagepath,experimentpath,'Raw\'];
datadir = [projectpath,experimentpath,'Data_bgtemplateoffsetcalc\'];
filetag={'DAPI','FITC','Cy3','Cy5'};
numrows=numel(rowmat);
numcols=numel(colmat);
numwells=numrows*numcols;
tempshot=[num2str(rowmat(1)),'_',num2str(colmat(1)),'_',num2str(sitemat(1))];
tempimg=single(imread([rawdir,tempshot,'_',filetag{1},'_stain.tif']));
for i=1:length(filetag)
    meanrawrel=ones(numel(tempimg),numel(sitemat))*NaN;
    for s=1:numel(sitemat)
        site=sitemat(s);
        highflierflag=1;
        hidx=0;
        while highflierflag
            hidx=hidx+1; rowidx=ceil(hidx/numcols); colidx=mod(hidx,numcols); colidx=colidx+(colidx==0)*numcols;
            refshot=[num2str(rowmat(rowidx)),'_',num2str(colmat(colidx)),'_',num2str(site)];
            refraw=single(imread([rawdir,refshot,'_',filetag{i},'_stain.tif']));
            [highflierflag,highthresh,refraw]=checkhighfliers(refraw,5);
        end
        %refraw=removesmears(refraw,highthresh);
        refrel=refraw(:)-nanmean(refraw(:));
        rawrel=ones(numel(refraw),numwells)*NaN;
        NCS=ones(numwells,1)*NaN;        
        cc=0;
        for row=rowmat
            for col=colmat
                cc=cc+1;
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [rawrel(:,cc),NCS(cc)]=getoffsetdiff(refrel,rawdir,shot,filetag{i},highthresh);
            end
        end
        goodwells=NCS>0.90;
        if sum(goodwells)<=2
            fprintf('Less than 3 good wells for file %s refshot %s\n',filetag{i},refshot);
            fprintf('sum(goodwells) = %0.0f\n',sum(goodwells));
        end
        rawrel=rawrel(:,goodwells);
        meanrawrel(:,s)=nanmean(rawrel,2);
    end
    filename=[datadir,'Background_',filetag{i},'.mat'];
    save(filename,'meanrawrel');
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