function [ mask ] = preprocessImage( imagePath )
    % Read the image
    image = imread(imagePath);

    % Convert data type
    I = double(image);

    %  0 - 1 scale
    I = mat2gray(I);

    % Apply adaptive thresholding
    mask = adaptivethresh(I, 10, -10, 'gaussian', 'relative');

    % Apply further filtering to clean salt-pepper noise & smaller objects
    mask = imfill(mask, 'holes');
    mask = imopen(mask, strel('disk', 4,0));
    mask = bwareaopen(mask, 100,4);
end
