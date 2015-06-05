function Average_ShadingImages
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadingpath='E:\20131126 R-point SR\20140510 46Hys-DoseResponse\Second Imaging\Serum\Raw\';

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names={
    %'BG';
    'DAPI';
    %'CFP';
	'YFP';
	'TxRed';
	'Cy5';
};
rowmat=4; %3 6
colmat=3:4; %3 9
sitemat=1:6;
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
    scarrays=[];
    for r=1:numrows
        row=rowmat(r);
        for c=1:numcols
            col=colmat(c);
            for s=1:numsites
                site=sitemat(s);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                block=single(imread([shadingpath,shot,'_',name,'_serum.tif']));
                if subtractcameranoise
                    block=block-bgcmos;
                end
                %%% For calculating normalization
                sctemp=shadingcorrection_1([],block,blurrad);
                %%% No normalization (to calculate camera noise)
                %sctemp=imfilter(block,fspecial('disk',blurrad),'symmetric');
                
                scarrays=cat(3,scarrays,sctemp);
                shadingcorrection=mean(scarrays,3);
                save([shadingpath,names{i},'_',num2str(site),'.mat'],'shadingcorrection');
            end
        end
    end
    %shadingcorrection=mean(scarrays,3);
    %save([shadingpath,names{i},'_',num2str(sitemat),'.mat'],'shadingcorrection');
    %save([shadingpath,names{i},'.mat'],'shadingcorrection');
end