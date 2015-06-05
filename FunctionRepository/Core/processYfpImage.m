function [ yfpprops ] = processYfpImage( inImage, mask )
    % Read the input image
    yfpImage = double(imread(inImage));
    yfpImageCorrected = illum_channel_block(yfpImage, imdilate(mask, strel('disk',5)));

    % Return properties of the objects
    yfpprops = regionprops(bwlabel(mask),yfpImageCorrected,'all');
end

