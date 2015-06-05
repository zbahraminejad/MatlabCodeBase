workspace;
close all;
clc;

parentpath = '/Users/Zahrabahrami/Desktop/Jan10th1215--Went pAKT_Plate_881-2/TimePoint_1';
parentList = dir(parentpath);parentList = {parentList.name};
allList = regexp(parentList, '.*[0-9].TIF', 'match');
allList = allList(~logical(cellfun(@isempty, allList)));
allList = cellfun(@cell2mat, allList,'UniformOutput',false);
channel = {'Dapi','YFP','mCherry'};

for rowi = 'A':'H'
    for coli = 4:4
        wellhit = regexp(allList, [rowi, sprintf('%0.2d',coli)]);
        wellPaths = allList(~logical(cellfun(@isempty, wellhit)));
    end
end

mask = preprocessImage('/Users/Zahrabahrami/Desktop/Jan10th1215--Went pAKT_Plate_881-2/TimePoint_1/Jan10th1215--Went pAKT_A04_w1.TIF');
imagesc(mask)
% convert binary image to numbered
imgPath = '/Users/Zahrabahrami/Desktop/Jan10th1215--Went pAKT_Plate_881-2/TimePoint_1/Jan10th1215--Went pAKT_A04_w2.TIF';
yfpprops = processYfpImage(imgPath, mask);
hist([yfpprops.MeanIntensity], 100)

%% Get a cytoplasmic mask
% The easy way is to use imdilate or bwmorph(bw, 'thicken')
% Use the nuclear mask as anti-mask. e.g. ~nuclearMask
dilatedBw= imdilate(label, strel('disk',4));
cytoBw = dilatedBw.*~label;

%% plot for example
hist([yfpprops.MeanIntensity], 100)