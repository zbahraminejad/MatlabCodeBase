function real=bgsubmin(raw,blurrad)
blur=imfilter(raw,fspecial('disk',blurrad),'symmetric');
c=min(blur(:));
bg=ones(size(raw))*c;
real=raw-bg;