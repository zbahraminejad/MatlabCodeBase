function [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_1(rawdir,name1,frames,nucr,debrisarea,badframes,RGBwrite,RGBdir)
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    raw1=log(single(imread([rawdir,name1,num2str(frames(i)),'.tif'])));
    %nuc_mask=blobdetector(raw1,nucr,blobthreshold,debrisarea);
    nuc_mask=threshmask(raw1,3);
    nuc_mask=markershed(nuc_mask,round(nucr*2/3));
    nuc_mask=bwareaopen(nuc_mask,debrisarea);
    nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
    
    [~,numcells(i)]=bwlabel(nuc_mask);
    nuc_area(i)=median(cell2mat(struct2cell(regionprops(nuc_mask,'Area'))));
end
dims=size(nuc_mask);
height=dims(1); width=dims(2);
blurthreshhigh=1.10*nanmedian(nuc_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=find(numcells>0 & nuc_area<blurthreshhigh,1,'first');
if firstgoodindex>1
    for i=1:firstgoodindex-1
        badframes(frames(1)-1+i)=1;
        if RGBwrite
            imwrite(uint16(zeros(dims)),[RGBdir,namenucedge,num2str(frames(i)),'.tif']);
        end
    end
end
badframes(frames(1)-1+firstgoodindex)=0;
blurthreshhigh=1.08*nuc_area(firstgoodindex);
blurthreshlow=0.95*nuc_area(firstgoodindex);
numthresh=0.8*numcells(firstgoodindex); %was 0.5
end

%{
%%% debugging: view images %%%%%%%%%%
[height,width]=size(raw1);
extractmask=bwmorph(nuc_mask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=extractmask;
tempframe(:,:,2)=imadjust(mat2gray(raw1));
tempframe(:,:,3)=nuc_mask;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(nuc_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
%}