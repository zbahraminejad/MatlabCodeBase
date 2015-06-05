function Immunostain_2(rowmat,colmat,sitemat)
fun=@(block_struct) nanmean(block_struct.data(:));  %function handle for blockproc call
a=single(imread('H:\Images\2013-12-12_BackgroundCharacterization\Raw\1_1_1_Cy3_stain.tif'));
a=imfilter(a,fspecial('disk',10));
[height,width]=size(a);
blocknum=11;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
%a=blockproc(a,[blockheight blockwidth],fun);
b=single(imread('H:\Images\2013-12-12_BackgroundCharacterization\Raw\2_10_1_Cy3_stain.tif'));
b=imfilter(b,fspecial('disk',10));
%b=blockproc(b,[blockheight blockwidth],fun);
% aorg=a;
% a=aorg(2:11,2:11);
% b=aorg(1:10,1:10);

%as=a(:); meana=mean(as); stda=std(as); norma=(as-meana)/stda; bs=b(:); meanb=mean(bs); stdb=std(bs); normb=(bs-meanb)/stdb; axb=norma.*normb; c=sum(axb(:))/(numel(as));
as=a(:)-mean(a(:)); bs=b(:)-mean(b(:)); c=dot(as/norm(as),bs/norm(bs));
fprintf('c = %8.4f\n',c);

imagesc(b); colorbar;
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%%% nucleoli
tempimage=mat2gray(YFP_raw); tempframe=imadjust(tempimage,stretchlim(tempimage,[0.01 0.9999]));
%}