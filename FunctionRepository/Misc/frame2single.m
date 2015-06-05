
imagepath = 'G:\Michael\';
experimentpath='20150401-CC-Diff\';


singledir = [imagepath,experimentpath,'Raw - Single\'];

rows = 2:7;
cols = [1,3,4,6];
sites = 1:2;
numrows = numel(rows);
numcols = numel(cols);
numsites = numel(sites);
channels = {'YFP','Cy5'};
numchannels = numel(channels);

frames = [113,245,377,509,653];
numframes = numel(frames);

for row = 1:numrows
    for col = 1:numcols
        for site = 1:numsites
            for frame = 1:numframes
                shot = [num2str(rows(row)),'_',num2str(cols(col)),'_',num2str(sites(site))];
                wellname = nameandsite(shot);
                rawdir=[imagepath,experimentpath,'Raw\',wellname];
                filename1 = [shot,'_',channels{1},'_',num2str(frames(frame))];
                filename2 = [shot,'_',channels{2},'_',num2str(frames(frame))];
                copyfile([rawdir,filename1,'.tif'],[singledir,filename1,'.tif'],'f');
                copyfile([rawdir,filename2,'.tif'],[singledir,filename2,'.tif'],'f');
            end
        end
    end
end