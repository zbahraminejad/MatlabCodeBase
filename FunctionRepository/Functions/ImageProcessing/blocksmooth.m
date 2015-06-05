function finalimg=blocksmooth(initimg,blocknum)
fun=@(block_struct) nanmedian(block_struct.data(:));
funnum=@(block_struct) sum(~isnan(block_struct.data(:)));
[height,width]=size(initimg);
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
blockimg=blockproc(initimg,[blockheight blockwidth],fun);
blocknumimg=blockproc(initimg,[blockheight blockwidth],funnum);
lowsample=blocknumimg<30;
totalmean=nanmean(blockimg(:));
if isnan(totalmean)
    fprintf('no open spaces!\n');
end
blockimg(lowsample)=totalmean;
finalimg=imresize(blockimg,[height width],'bicubic');
