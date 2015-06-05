projectpath='H:\Documents\Projects\';
imagepath='H:\Images\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
datadir=([projectpath,experimentpath,'Data\']);
rawdir=[imagepath,experimentpath,'Raw\'];
nucname='CFP_';
primaryname='TexasRed_';
secondaryname='Cy5_';
nucr=12;

%rows={'B'};cols={'03'};sites={'1'};
%%% 20130715 wells %%%%%%%%%%%%%%%
rows={'B' 'C' 'D' 'E'};
cols={'03' '04' '05' '06'};
sites={'1' '2' '3' '4'};
%%% 20130719 wells %%%%%%%%%%%%%%%
%rows={'B' 'C' 'D' 'E'};
%cols={'03' '11'};
%sites={'1' '2' '3' '4'};
%%% 20130905 wells %%%%%%%%%%%%%%%
%rows={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
%cols={'02'};
%sites={'1' '2' '3' '4'};
stain='prestain';

totalwells=length(rows)*length(cols)*length(sites);
bleedthroughrates=zeros(totalwells,2);
indexcount=0;
for rowidx=1:length(rows)
    row=rows{rowidx};
    for colidx=1:length(cols)
        col=cols{colidx};
        for siteidx=1:length(sites)
            indexcount=indexcount+1;
            site=sites{siteidx};
            shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
            
            nucfile=[rawdir,shot,'\',nucname,stain,'.tif'];
            primaryfile=[rawdir,shot,'\',primaryname,stain,'.tif'];
            secondaryfile=[rawdir,shot,'\',secondaryname,stain,'.tif'];
            nuc_raw=single(imread(nucfile));
            primary_raw=single(imread(primaryfile));
            secondary_raw=single(imread(secondaryfile));
            nuc_mask=blobdetector(log(nuc_raw),nucr,-0.02);
            nuc_center=squeeze(cell2mat(struct2cell(regionprops(nuc_mask,'Centroid')')))';
            [primary_bg,secondary_bg]=getbackground_12(nuc_mask,nuc_center,primary_raw,secondary_raw,nucr,10,0.25,50,50);
            primary_raw=cell2mat(struct2cell(regionprops(nuc_mask,primary_raw,'MeanIntensity')))';
            secondary_raw=cell2mat(struct2cell(regionprops(nuc_mask,secondary_raw,'MeanIntensity')))';
            primary_bs=primary_raw-primary_bg;
            secondary_bs=secondary_raw-secondary_bg;
            primary_bs=primary_bs(:);
            secondary_bs=secondary_bs(:);
            p=robustfit(primary_bs,secondary_bs);
            %dscatter(primary_bs,secondary_bs);
            %xmin=min(primary_bs); xmax=max(primary_bs); xstep=(xmax-xmin)/100; xvals=xmin:xstep:xmax;
            %hold on; plot(xvals,p(1)+p(2)*xvals,'r');
            bleedthroughrates(indexcount,:)=[p(1) p(2)]; %store [intercept slope]
        end
    end
end
bleedthroughrate=mean(bleedthroughrates);
save([datadir,'bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');