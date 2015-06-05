function finalimg=blocksmooth(initimg,blocknum)
fun=@(block_struct) nanmedian(block_struct.data(:));
[height,width]=size(initimg);
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
blockimg=blockproc(initimg,[blockheight blockwidth],fun);
finalimg=imresize(blockimg,[height width],'bicubic');
