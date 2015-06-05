function bgsubbed=bgsub_MC_mask(org,mask,blockheight,blockwidth)
fun=@(block_struct) nanmedian(block_struct.data(:));  %function handle for blockproc call

bgorg=org.*~mask;
bgorg(bgorg==0)=NaN;
bgblock=blockproc(bgorg,[blockheight blockwidth],fun);
bgfinal=imresize(bgblock,size(org),'bicubic');
bgsubbed=org-bgfinal;
%bgsubbed(bgsubbed<0)=0;  %previously was commented out