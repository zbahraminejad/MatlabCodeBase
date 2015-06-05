function finalimg=blocksmooth_mode_3(initimg,blocknum,globalbg)
fun=@(block_struct) calcmode(block_struct.data(:),globalbg);
[height,width]=size(initimg);
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
blockimg=blockproc(initimg,[blockheight blockwidth],fun);
finalimg=imresize(blockimg,[height width],'bicubic');
end

function modeval=calcmode(vals,globalmedian)
if sum(~isnan(vals))<100
    modeval=globalmedian;
else
    binmax=prctile(vals,95);
    binmin=prctile(vals,5);
    vals=vals(vals<binmax & vals>binmin);
    [kval,xval]=ksdensity(vals);
    modeval=xval(kval==max(kval));
    modeval=modeval(1);
end
end