function RGB=randomlabel(labelimg,numcolors)
colormapsize=numcolors;
colormap=jet(colormapsize); close(gcf);
numobs=max(labelimg(:));
colormapidx=ceil(rand(numobs,1)*colormapsize);
cellcolor=colormap(colormapidx,:);
[height,width]=size(labelimg);
emptyframe=zeros(height,width);
R=emptyframe; G=emptyframe; B=emptyframe;
for i=1:numobs
    labelidx=labelimg==i;
    R(labelidx(:))=cellcolor(i,1);
    G(labelidx(:))=cellcolor(i,2);
    B(labelidx(:))=cellcolor(i,3);
end
RGB=zeros(height,width,3);
RGB(:,:,1)=R;
RGB(:,:,2)=G;
RGB(:,:,3)=B;
RGB(RGB>1)=1;
end