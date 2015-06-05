function BackgroundSurvey_CrossCorr
rowmat=1:7;colmat=1:12;sitemat=1;
%rowmat=1;colmat=1:2:11;sitemat=1;
refshot='1_1_1';
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
experimentpath='2013-12-12_BackgroundCharacterization\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
rawdir = [imagepath,experimentpath,'Raw\'];
filetag='Cy3_';

refraw=single(imread([rawdir,refshot,'_',filetag,'stain.tif']));
highthresh=callthresh(refraw);
refraw=removesmears(refraw,highthresh);
refrel=refraw(:)-nanmean(refraw(:));

NCS=ones(numel(rowmat),numel(colmat),numel(sitemat))*NaN; offsetdiff=NCS;
cc=0;
for row=rowmat
    for col=colmat
        for site=sitemat
            cc=cc+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            [NCS(row,col,site),offsetdiff(row,col,site)]=calcCC(refrel,rawdir,shot,filetag,highthresh);
        end
    end
end
%imagesc(colmat,rowmat,NCS(:,1:12,4),[0.75 1.00]); colorbar; set(gca,'XTick',[1:1:12]);
imagesc(colmat,rowmat,offsetdiff(:,1:12,1),[0 15]); colorbar; set(gca,'XTick',[1:1:12]);
set(gcf,'color','w','PaperPosition',[0 0 8 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end

function [NCS,offsetdiff]=calcCC(refrel,rawdir,shot,filetag,highthresh)
filename=[rawdir,shot,'_',filetag,'stain.tif'];
raw=single(imread(filename));
raw=removesmears(raw,highthresh);
rawrel=raw(:)-nanmean(raw(:));
nanidx=isnan(refrel) | isnan(rawrel);
refrel(nanidx)=[]; rawrel(nanidx)=[];
NCS=dot(refrel/norm(refrel),rawrel/norm(rawrel));
offsetdiff=mean(abs(refrel-rawrel));
end

%{
%%% debugging: view images %%%%%%%%%%
imagesc(rawimage); colorbar;
figure; imshow(mask);
%}