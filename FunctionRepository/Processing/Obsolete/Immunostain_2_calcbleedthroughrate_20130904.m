projectpath = 'H:\Documents\Projects\';
imagepath = 'H:\Images\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
datadir = ([projectpath,experimentpath,'Data\']);
rawdir = [imagepath,experimentpath,'Raw\'];
RFPname = 'TexasRed_'; %p21
IFname = 'Cy5_';
nucr=12;

%rows={'B'};cols={'03'};sites={'1'};
%%% 20130715 wells %%%%%%%%%%%%%%%
rows = {'B' 'C' 'D' 'E'};
cols = {'03' '04' '05' '06'};
sites = {'1' '2' '3' '4'};
%%% 20130719 wells %%%%%%%%%%%%%%%
% rows = {'B' 'C' 'D' 'E'};
% cols = {'03' '11'};
% sites = {'1' '2' '3' '4'};

totalwells=length(rows)*length(cols)*length(sites);
bleedthroughrates=zeros(totalwells,2);
indexcount=0;
for rowidx = 1:length(rows)
    row = rows{rowidx};
    for colidx = 1:length(cols)
        col = cols{colidx};
        for siteidx = 1:length(sites)
            indexcount=indexcount+1;
            site = sites{siteidx};
            shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
            RFPfile=[rawdir,shot,'\',RFPname,'prestain.tif'];
            IFfile=[rawdir,shot,'\',IFname,'prestain.tif'];
            RFP_raw=single(imread(RFPfile));
            IF_raw=single(imread(IFfile));
            RFP_bs=bgsub(log(RFP_raw),10*nucr,0.05);
            IF_bs=bgsub(log(IF_raw),10*nucr,0.05);
            RFP_bs=RFP_bs(:);
            IF_bs=IF_bs(:);
            
            %outeridx=find(RFP_bs<prctile(RFP_bs,1) | RFP_bs>prctile(RFP_bs,99));
            %RFP_outer=RFP_bs(outeridx); IF_outer=IF_bs(outeridx);
            p=robustfit(RFP_bs,IF_bs);
            %p=fitline(RFP_outer,IF_outer);
            %p=fitlineTotalLeastSquares(RFP_bs(:),IF_bs(:));
            %scatter(RFP_outer,IF_outer);
            %xvals=[min(RFP_outer);max(RFP_outer)];yvals=xvals*p.m+p.b;
            %hold on; plot(xvals,yvals,'r');
            density_scatter_heatmap(RFP_bs,IF_bs,200,200);
            
            bleedthroughrates(indexcount,:)=[p.m p.b]; %store slope and intercept
        end
    end
end
bleedthroughrate=mean(bleedthroughrates);
save([datadir,'bleedthroughrate.mat'],'bleedthroughrate');