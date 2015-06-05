function [firstgoodindex,blurthresh,numthresh,badframes,height,width]=timelapsesetup(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite)
%%% determine median cell size for blur detection %%%%%%%%%%%%%%%%%%%%%%%%%
nuc_area=zeros(3,1); numcells=zeros(3,1);
for i=1:3
    raw1=log(single(imread([rawdir,name1,num2str(frames(i)),'.tif'])));
    nuc_mask=blobdetector(raw1,nucr,blobthreshold,debrisarea);
    [~,numcells(i)]=bwlabel(nuc_mask);
    nuc_area(i)=median(cell2mat(struct2cell(regionprops(nuc_mask,'Area'))));
end
dims=size(nuc_mask);
height=dims(1); width=dims(2);
blurthresh=1.10*median(nuc_area);
%%% determine first good frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstgoodindex=find(numcells>0 & nuc_area<blurthresh,1,'first');
if firstgoodindex>1
    for i=1:firstgoodindex-1
        badframes(i)=1;
        if maskwrite
            imwrite(uint16(zeros(dims)),[maskdir,namenucedge,num2str(frames(i)),'.tif']);
        end
    end
end
badframes(firstgoodindex)=0;
blurthresh=1.08*nuc_area(firstgoodindex);
numthresh=0.5*numcells(firstgoodindex);
end