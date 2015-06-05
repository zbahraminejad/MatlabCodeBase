function CompressImages(row,col,site)

imagepath='H:\Images\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
orgdir=[imagepath,experimentpath,'Raw\',shot,'\'];
compdir=[imagepath,experimentpath,'Raw_Compressed\',shot];
if ~exist(compdir,'dir')
    mkdir(compdir);
end

SF=1;EF=230;
nucname='CFP_'; %nuc
YFPname='YFP_'; %DHB
RFPname='TexasRed_'; %p21

timetotal=tic;
%%% derive cell trace signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f=SF:EF
    timeframe=tic;
    nucfile=[nucname,num2str(f),'.tif'];
    YFPfile=[YFPname,num2str(f),'.tif'];
    RFPfile=[RFPname,num2str(f),'.tif'];
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_raw=single(imread([orgdir,nucfile]));
    YFP_raw=single(imread([orgdir,YFPfile]));
    RFP_raw=single(imread([orgdir,RFPfile]));
    %%% compress images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_comp=imresize(nuc_raw,0.5,'bilinear');
    YFP_comp=imresize(YFP_raw,0.5,'bilinear');
    RFP_comp=imresize(RFP_raw,0.5,'bilinear');
    %%% save images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imwrite(uint16(nuc_comp),[compdir,'\',nucfile]);
    imwrite(uint16(YFP_comp),[compdir,'\',YFPfile]);
    imwrite(uint16(RFP_comp),[compdir,'\',RFPfile]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
toc(timetotal);