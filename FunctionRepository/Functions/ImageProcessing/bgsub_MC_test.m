function bgfinal=bgsub_MC(org,blockheight,blockwidth)
fun=@(block_struct) ThreshImage_MC_block(block_struct.data(:));  %function handle for blockproc call

bgblock=blockproc(org,[blockheight blockwidth],fun);
bgfinal=imresize(bgblock,size(org),'bicubic');
bgsubbed=org-bgfinal;
%bgsubbed(bgsubbed<0)=0;  %previously was commented out