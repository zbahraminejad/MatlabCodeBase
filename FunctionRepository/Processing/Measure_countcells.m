%imagepath='H:\Images\';
imagepath='E:\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
experimentpath='20131213 R-point CC\20140423 3C G1S 46i\';
nucr=12;
nucname='H2B_';

rows={'2' '3' '4' '5' '6' '7'};
cols={'2' '3' '4' '5' '6' '7' '8' '9' '10' '11'};
sites={'1' '2' '3' '4'};
frame=1;

numsites=numel(rows)*numel(cols)*numel(sites);
cellcount=cell(numsites,2);
i=0;
for rowidx=1:length(rows)
    row=rows{rowidx};
    for colidx=1:length(cols)
        col=cols{colidx};
        for siteidx=1:length(sites)
            i=i+1;
            site=sites{siteidx};
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            cellcount{i,1}=shot;
            rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
            raw1=single(imread([rawdir,nucname,num2str(frame),'.tif']));
            [nuc_mask,~]=blobdetector_foreground(log(raw1),nucr,-0.02,200);
            [~,numcells]=bwlabel(nuc_mask);
            cellcount{i,2}=numcells;
        end
    end
end
keyboard;