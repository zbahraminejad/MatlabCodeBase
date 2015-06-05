function markcellsbygate(datadir,dapilist,dapidir,centroids,cebpbvalues,ppargvalues,cebpbGate,ppargGate,markradius,signum)

dummy = double(imread([dapidir,'\',dapilist{1}]));
[height,width] = size(dummy);
emptyframe = zeros(height,width);
R = emptyframe;G=emptyframe;B=emptyframe;

for i = 1:length(dapilist)

dapi = double(imread([dapidir,'\',dapilist{i}]));
dapi = imadjust(mat2gray(dapi));
cents = round(centroids{i});
data = [cents ppargvalues{i}' cebpbvalues{i}'];
testx = data(:,1);
testy = data(:,2);
centroidindex = sub2ind([height,width],testy,testx);
centroidmask = emptyframe;
centroidmask(centroidindex) = 1;
centroidmask = imdilate(centroidmask,strel('disk',markradius,0));
dapi(centroidmask==1) = 0;

highcentroids = data(data(:,4)> cebpbGate & data(:,3)>ppargGate,1:2);
highind = sub2ind([height,width],round(highcentroids(:,2)),round(highcentroids(:,1)));
highmask = emptyframe;
highmask(highind)=1;

R = imdilate(highmask,strel('disk',markradius,0));

lowcentroids = data(data(:,4)< cebpbGate & data(:,3)<ppargGate,1:2);
lowind = sub2ind([height,width],lowcentroids(:,2),lowcentroids(:,1));
lowmask = emptyframe;
lowmask(lowind) = 1;

G = imdilate(lowmask,strel('disk',markradius,0));

weirdcentroids = data(data(:,4)> cebpbGate & data(:,3)<ppargGate,1:2);
weirdind = sub2ind([height,width],weirdcentroids(:,2),weirdcentroids(:,1));
weirdmask = emptyframe;
weirdmask(weirdind) = 1;

B = imdilate(weirdmask,strel('disk',markradius,0));

outcentroids = data(data(:,4)< cebpbGate & data(:,3)>ppargGate,1:2);
outind = sub2ind([height,width],outcentroids(:,2),outcentroids(:,1));
outmask = emptyframe;
outmask(outind) = 0.75;

R2 = imdilate(outmask,strel('disk',markradius,0));
G2 = imdilate(outmask,strel('disk',markradius,0));
B2 = imdilate(outmask,strel('disk',markradius,0));

finalimage = zeros(height,width,3);
finalimage(:,:,1) =   dapi + R + R2;
finalimage(:,:,2) =   dapi + G + G2;
finalimage(:,:,3) =   dapi + B + B2;
% figure('Name',dapilist{i}(1:end-9)),
% imshow(finalimage)

[classified,finaldistance,finalindex,highnum,lownum,weirdnum]=classifycell(cents,centroidindex,highind,lowind,weirdind,signum);
save([datadir,'classified\',dapilist{i}(1:end-7),'_','classified.mat'],'classified','finaldistance','finalindex','finalimage','highnum','lownum','weirdnum','finalimage');
end

end