function Average_ShadingImages
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath='D:\Images\';
shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20130607 p21dCy2\20140629 DHB_p21dCy1dK_SRvsCC\';
separatedirectories=1;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names={
    'H2B';
    'DHB';
	%'p21dCy1dK';
	%'pRb';
	%'Cy5';
};
rowmat=[4 7]; %[3 4 7]
colmat=[5:7]; %[5:7]
sitemat=1:3; %1
framemat=[1];
nucr=12;
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat); numframes=numel(framemat);
load([shadingpath,'BG.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
for i=1:length(names)
    name=char(names{i});
    for s=1:numsites
        rawarrays=[];
        site=sitemat(s);
        for c=1:numcols
            col=colmat(c);
            for r=1:numrows
                row=rowmat(r);
                for f=1:numframes
                    frame=framemat(f);
                    wellcolname=wellcol(row,col);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    if separatedirectories
                        rawdir=[imagepath,experimentpath,'Raw\'];
                        %rawdir=[imagepath,experimentpath,'Raw\',wellcolname,'\site_',num2str(site),'\'];
                        %raw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        raw=single(imread([rawdir,shot,'\',names{i},'_',num2str(frame),'.tif']));
                    else
                        rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
                        raw=single(imread([rawdir,names{i},'_stain.tif']));
                    end
                    raw=raw-bgcmos;
                    raw=imfilter(raw,fspecial('disk',5),'symmetric');
                    rawarrays=cat(3,rawarrays,raw);
                end
            end
        end
        rawcalc=min(rawarrays,[],3);
        %rawcalc=mean(rawarrays,3);
        inferredbg=shadingcorrection_1([],rawcalc,10*nucr); %can be 1*nucr with enough sites
        save([imagepath,experimentpath,'Raw\IBG\',names{i},'_',num2str(site),'.mat'],'inferredbg');
        %save([imagepath,experimentpath,'Raw\IBG_',names{i},'.mat'],'inferredbg');
    end
end