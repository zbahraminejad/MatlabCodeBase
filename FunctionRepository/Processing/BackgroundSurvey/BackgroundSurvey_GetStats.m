function BackgroundSurvey_getstats
rowmat=1:7;colmat=1:12;sitemat=1:4;
%rowmat=1;colmat=3;sitemat=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath = 'H:\Images\';
experimentpath='2013-12-12_BackgroundCharacterization\';
%experimentpath='2013-11-26_CycDp21SerumRelease\Experiment_20131216\';
filetag='Cy3_';
rawdir = [imagepath,experimentpath,'Raw\'];
refshot=[num2str(rowmat(1)),'_',num2str(colmat(1)),'_',num2str(sitemat(1))];
refimg=single(imread([rawdir,refshot,'_',filetag,'stain.tif']));
highthresh=callthresh(refimg);
meanint=ones(numel(rowmat),numel(colmat),numel(sitemat))*NaN; stdint=meanint;
cc=0;
for row=rowmat
    for col=colmat
        for site=sitemat
            cc=cc+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            [meanint(row,col,site),stdint(row,col,site)]=getstats(rawdir,shot,filetag,highthresh);
        end
    end
end
%scatter(meanint(:),stdint(:));
%xlabel('mean RFU'); ylabel('stddev RFU');
imagesc(1:12,1:7,stdint(:,:,1)); colorbar;
set(gcf,'color','w','PaperPosition',[0 0 8 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end

function [meanint,stdint]=getstats(rawdir,shot,filetag,highthresh)
filename=[rawdir,shot,'_',filetag,'stain.tif'];
raw=single(imread(filename));
raw=removesmears(raw,highthresh);
meanint=nanmean(raw(:));
stdint=nanstd(raw(:));
end