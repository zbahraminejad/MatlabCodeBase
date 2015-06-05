function maskwatershed=markershed(mask,raw,nucrad)
valleys=-bwdist(~mask);
blur=imfilter(raw,fspecial('disk',round(nucrad/4)),'symmetric');
blur=imfilter(blur,fspecial('gaussian',nucrad*2,nucrad),'symmetric');
basins=imregionalmax(blur);
basins=imdilate(basins,strel('disk',round(nucrad/4)));
bigvalleys=bwdist(mask);
outerridges=watershed(bigvalleys);
outerridges=outerridges==0;
finalvalleys=imimposemin(valleys,basins | outerridges);
finalridges=watershed(finalvalleys);
mask(finalridges==0)=0;
maskwatershed=mask;
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(mask,'remove');
[height,width]=size(mask);
RGB=zeros(height,width,3);
RGB(:,:,1)=basins;
RGB(:,:,2)=largecellmarker;
RGB(:,:,3)=extractmask;
figure,imshow(RGB);
%}