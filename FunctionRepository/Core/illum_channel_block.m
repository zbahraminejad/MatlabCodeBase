% 5x5 block mode filter. You can give a binary mask if you want to exclude
% that region from the caluclation of mode.

function image = illum_channel_block(image, mask)

    blocknum = 5;
    
    [height, width] = size(image);
    blockheight = ceil(height/blocknum);
    blockwidth = ceil(width/blocknum);
    
    tempImage = image;
    if nargin == 2
        tempImage(mask) = NaN;
    end
    fun = @(block_struct) blockfun_ksdensity(block_struct);
    blockimg = blockproc(tempImage,[blockheight blockwidth],fun);
    background = imresize(blockimg,[height width],'bicubic');
    %%
    image = image - background;
    image = image - min(image(:));

end

function [modePixel] = blockfun_ksdensity(block_struct)
    pixels = block_struct.data(:);
    pixels = pixels(pixels<prctile(pixels,95) & pixels>prctile(pixels,5));
    [kval,xval] = ksdensity(double(pixels));
    modePixel = xval(kval == max(kval));
end

