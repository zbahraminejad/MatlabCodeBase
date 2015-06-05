function bg=getbg_simple(img,rd,bgpeak)
% bgpeak: fraction (0-1); if zero, will be generated automatically

%% calculate the ratio
ratio=min([1,20/rd]);
rd_new=min([rd,20]);
mkbg=getnhood(strel('disk',rd_new,0));szb=size(mkbg,1);
mksm=getnhood(strel('disk',round(rd_new/2),0));szs=size(mksm,1);
x1=1+(szb-szs)/2;x2=x1+szs-1;
mk=mkbg;mk(x1:x2,x1:x2)=mkbg(x1:x2,x1:x2)-mksm;

%% calculate bg peak
if bgpeak==0
    bgpeak=(4-pi)/8;
end

%% use ordfilt2 to get bg
im_small=imresize(img,ratio,'bicubic');
bg_small=ordfilt2(im_small,round(sum(mk(:))*bgpeak),mk,'symmetric');
bg0=imresize(bg_small,size(img),'bicubic');
bg=imfilter(bg0,fspecial('disk',min([rd,32])),'symmetric');