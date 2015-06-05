function bias=illumbiasEmptyWell(raw,noise)
% Blur image
blur=imfilter(raw,fspecial('disk',10),'symmetric'); % change disk size back to 5
% Subtract camera noise from image
backgroundonly=blur-noise;
% % Remove foreground pixels by setting them to NaN (Not a Number)
% nanmask=imdilate(mask,strel('disk',5));
% backgroundonly(nanmask)=NaN;
% Break the image into sections and estimate the background in each section
% by a given intensity percentile. Interpolate all of the sections to
% return a smoothened background image.
blocknum=5;
prctilethresh=10;
background=blockpercentile(backgroundonly,blocknum,prctilethresh);
maxillumination=max(background(:));
bias=background/maxillumination;
end

function prctileimage=blockpercentile(initimg,blocknum,prctilethresh)
getpercentile=@(block_struct) prctile(block_struct.data(:),prctilethresh);
[height,width]=size(initimg);
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
blockimg=blockproc(initimg,[blockheight blockwidth],getpercentile);
prctileimage=imresize(blockimg,[height width],'bicubic');
end