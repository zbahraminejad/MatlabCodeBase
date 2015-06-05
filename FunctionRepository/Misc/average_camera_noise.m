noisedir = 'D:\Michael\Camera Noise_1080_bin2\';
shadingpath='D:\Michael\ShadingImages\DCYTC 10x\';
numreps = 5;
resolution = 1080;
filetitle = 'CameraNoise_1080_bin2';

storemat = zeros(resolution,resolution,numreps);

for rep = 1:numreps;
    bg = double(imread([noisedir,filetitle,'_',num2str(rep),'.tif']));
    storemat(:,:,rep) =  bg;
end

averagebg = mean(storemat,3);
BG = imfilter(averagebg,fspecial('average',3));

save([shadingpath,filetitle,'.mat'],'BG');