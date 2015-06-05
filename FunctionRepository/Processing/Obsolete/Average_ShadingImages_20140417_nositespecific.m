function Average_ShadingImages
%row=1;col=1;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shadingpath= 'H:\Images\ShadingImages\20140402_DAPI_CFP_YFP_TxRed_Cy5\';

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names={
    'DAPI_50ms_pink';
	'CFP_25ms_pink';
	'YFP_25ms_yellow';
	'TxRed_50ms_red';
	'Cy5_100ms_red';
};
rowmat=[3 6]; %3 6
colmat=[3 9]; %3 9
sitemat=1:4;
blurrad=10; %pixelradius of blurring filter
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat);
%sccell=cell(numrows,numcols,numsites);
shadingcorrection=cell(length(names));
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
                block=single(imread([shadingpath,'Raw\',shot,'_',name,'.tif']));
                %sccell{r,c,s}=shadingcorrection_1([],block,blurrad);
                sctemp=shadingcorrection_1([],block,blurrad);
                scarrays=cat(3,scarrays,sctemp);
            end
        end
    end
    shadingcorrection=mean(scarrays,3);
    save([shadingpath,names{i},'.mat'],'shadingcorrection');
end