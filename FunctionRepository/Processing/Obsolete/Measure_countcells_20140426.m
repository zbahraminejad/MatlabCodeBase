imagepath='H:\Images\';
experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
nucr=8;
nucname='CFP_';

rows={'B' 'C' 'D' 'E' 'F' 'G'};
cols={'02' '03' '04' '05' '06' '07' '08' '09' '10'};
sites={'1' '2'};

for rowidx=1:length(rows)
    row=rows{rowidx};
    for colidx=1:length(cols)
        col=cols{colidx};
        for siteidx=1:length(sites)
            site=sites{siteidx};
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
            nuc_raw=log(single(imread([rawdir,nucname,num2str(frame),'.tif'])));
            nuc_mask=blobdetector(nuc_raw,nucr,-0.03);
            [~,numcells]=bwlabel(nuc_mask);
        end
    end
end