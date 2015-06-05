projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
experimentpath='20131209 Pipdeg Optimization\20140402 PFACS\';
datadir=([projectpath,experimentpath,'Data\']);
rawdir=[imagepath,experimentpath,'Raw\'];
shadingpath= 'H:\Images\ShadingImages\20140402_DAPI_CFP_YFP_TxRed_Cy5\';
nucname='Hoechst_';
primaryname='p21dCy1_';
secondaryname='EdU_';
nucr=12; blobthreshold=-0.02; debrisarea=100;
separatedirectories=0;
if separatedirectories
    strvar='\';
else
    strvar='_';
end

%%% wells %%%%%%%%%%%%%%%%%%%%%%%%
rows={'1' '2' '3' '4' '5' '6' '7' '8'};
cols={'8'};
sites={'1' '2' '3' '4'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stain='stain';

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
            
            nucfile=[rawdir,shot,strvar,nucname,stain,'.tif'];
            primaryfile=[rawdir,shot,strvar,primaryname,stain,'.tif'];
            secondaryfile=[rawdir,shot,strvar,secondaryname,stain,'.tif'];
            nuc_raw=single(imread(nucfile)); %load([shadingpath,'DAPI.mat'],'shadingcorrection'); nuc_raw=nuc_raw./shadingcorrection;
            primary_raw=single(imread(primaryfile)); %load([shadingpath,'TxRed.mat'],'shadingcorrection'); primary_raw=primary_raw./shadingcorrection;
            secondary_raw=single(imread(secondaryfile)); %load([shadingpath,'Cy5.mat'],'shadingcorrection'); secondary_raw=secondary_raw./shadingcorrection;
            [nuc_mask,~]=blobdetector_foreground(log(nuc_raw),nucr,blobthreshold,debrisarea);
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