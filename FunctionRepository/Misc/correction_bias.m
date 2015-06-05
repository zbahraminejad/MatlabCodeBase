noisedir = 'D:\Michael\Camera Noise_10062014\';

filelist  = dir(noisedir);
filelist = {filelist.name};
filelist(1:2) = [];
dummy = imread([noisedir,filelist{1}]);
[height,width] = size(dummy);
reps = 5;
all = zeros(height,width,reps);
for i = 1:reps;
   all(:,:,i) = double(imread([noisedir,filelist{1}])) ;
end

meannoise = sum(all,3);
meannoise = meannoise./reps;

smoothnoise = imfilter(meannoise,fspecial('disk',3));

save('D:\Michael\ShadingImages\DCYTC 10x\cameranoise_4sites.mat','smoothnoise');