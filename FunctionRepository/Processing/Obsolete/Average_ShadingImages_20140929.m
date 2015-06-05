function Average_ShadingImages
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadingpath='E:\20140530-Henny-Plate2\Raw\Serum\';

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names={
    %'BG';
    %'DAPI';
    'YFP';
	'RFP';
	%'TxRed';
	'Cy5';
};
rowmat=3; %3 6
colmat=5:6; %3 9
sitemat=1:9;
blurrad=10; %pixelradius of blurring filter
%%% get camera noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subtractcameranoise=1;
if subtractcameranoise
    bgcmospath='H:\Images\ShadingImages\20140425 DCYTC 10x\';
    %bgcmospath='H:\Images\ShadingImages\20140522 20xBin2\';
    load([bgcmospath,'BG','.mat'],'shadingcorrection');
    bgcmos=shadingcorrection;
end
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat);
for i=1:length(names)
    name=char(names{i});
    for s=1:numsites
        scarrays=[];
        site=sitemat(s);
        for c=1:numcols
            col=colmat(c);
            for r=1:numrows
                row=rowmat(r);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                block=single(imread([shadingpath,shot,'_',name,'_stain.tif']));
                if subtractcameranoise
                    block=block-bgcmos;
                end
                %%% For calculating normalization
                sctemp=shadingcorrection_1([],block,blurrad);
                %%% No normalization (to calculate camera noise)
                %sctemp=imfilter(block,fspecial('disk',blurrad),'symmetric');
                
                scarrays=cat(3,scarrays,sctemp);
                
            end
        end
        shadingcorrection=mean(scarrays,3);
        save([shadingpath,names{i},'_',num2str(site),'.mat'],'shadingcorrection');
    end
end