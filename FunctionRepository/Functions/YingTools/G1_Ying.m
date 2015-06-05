% G1 tracing with better workflow
%% define path, make folder

path = pwd;

rawdir = [path,'/Raw/'];
datadir = [path,'/Data2/'];
mkdir(datadir);

%% input image names
V1 = 'H2B'; %V1 is used for nucleus identification and alignment
V2 = 'DHB'; %V2 for cyto/nucleus ratio
V3 = 'Geminin'; %V3 for intensity inside cell

lo = 0.5;
hi = 0.75;

%% load and initiate

V1bias = cell2mat(struct2cell(load(strcat(V1, '_bias.mat'),'bias')));
V2bias = cell2mat(struct2cell(load(strcat(V2, '_bias.mat'),'bias')));
V3bias = cell2mat(struct2cell(load(strcat(V3, '_bias.mat'),'bias')));

frames = 1:140;
totalframes = length(frames);
noise=double(imread([rawdir,'Noise.tif']));
maxcellnum=1000;
celldata=ones(maxcellnum,totalframes,6)*NaN; 
%%%centroid (x2), nucleus area, DHB nucleus intensity, DHB cyto intensity,
%%%Geminin intensity. intensity are calculated as the average between
%%%quantile lo and quantile hi.
alignment = zeros(totalframes,2); 

[height, width] = size(noise);

%% remove camera noise, flatten bias, remove background, align and extract for all

for f=frames
    fprintf('frame %0.0f\n',f);
    tic
    rawV1=double(imread([rawdir, V1,'_',num2str(f),'.tif']));
    
    %%% remove camera noise, correct for bias, minus exp background for V1  
    
    rawV1flat=(rawV1-noise)./V1bias;
    nucmaskV1=threshmask2(rawV1, 3);
    foreground=nucmaskV1;
    realV1=subtractbg(rawV1flat,foreground);
    
    %%% segment nuclei by watershed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nucmaskwatershedV1=markershed2(nucmaskV1, 8);
    nucmaskV1=bwareaopen(nucmaskwatershedV1,200);
    nucmaskV1=imclearborder(nucmaskV1,4);
    % assessmask(realV1, nucmaskV1);

    %%% remove camera noise, correct for bias, minus exp background for V2, V3  

    rawV2=double(imread([rawdir, V2, '_',num2str(f),'.tif']));
    rawV2flat=(rawV2-noise)./V2bias;
    nucmaskV2=threshmask2(rawV2, 3);
    realV2=subtractbg(rawV2flat,nucmaskV2);

    rawV3=double(imread([rawdir, V3, '_',num2str(f),'.tif']));
    rawV3flat=(rawV3-noise)./V3bias;
    nucmaskV3=threshmask(rawV3);
    realV3=subtractbg(rawV3flat,nucmaskV3);

    %%% get cytomask mask with the same label as nuclabel %%%%%%%%%%%%%%
    
    nuclabel = bwlabel(nucmaskV1);
    cellmask = bwmorph(nucmaskV1, 'thicken', 10);
    
    cellCC = bwconncomp(cellmask);
    numObjs = cellCC.NumObjects;
    imageSize = cellCC.ImageSize;
    celllabel = zeros(imageSize);
    for k = 1:numObjs
        IdxList = cell2mat(cellCC.PixelIdxList(k));
        [x, y] = ind2sub(imageSize, IdxList);
        xcenter = floor(mean(x));
        ycenter = floor(mean(y));
        celllabel(IdxList) = nuclabel(xcenter, ycenter);
    end
    % the following three lines are to make sure my labeling is correct
    % checkmask = celllabel .* nucmaskV1;
    % A = abs(checkmask - nuclabel);
    % sum(A(:)); sum(A(:))=0, correct label
    
    cytolabel = celllabel - nuclabel;
    foregrounddilate = imdilate(foreground, strel('disk', 2, 0));
    cytolabel(foregrounddilate == 1) = 0;
    % assessmask(realV2, cytomask);

    
    %%% Save nucleus mask and cyto mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imwrite(nucmaskV1, [datadir, 'nucmask_', num2str(f), '.tif']);
    imwrite(cytolabel, [datadir, 'cytomask_', num2str(f), '.tif']);

    % save logical mask is about 1/5 the time as saving a labelled mask %%%
    
    %%% Align %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if f>1
        %crosscorrscore=abs(normxcorr2(nucmaskprev,nucmask));
        crosscorrscore=abs(normxcorr2(foregroundprev,foreground));
        [~,idx]=max(crosscorrscore(:));
        [y,x]=ind2sub(size(crosscorrscore),idx);
        alignxrel=width-x;
        alignyrel=height-y;
        alignment(f,:)=alignment(f-1,:)+[alignxrel,alignyrel];
     end

    foregroundprev=foreground;

    
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nuc_info=struct2cell(regionprops(nucmaskV1,realV1,'Area','Centroid')');
    area=squeeze(cell2mat(nuc_info(1,1,:)));
    center=squeeze(cell2mat(nuc_info(2,1,:)))';
    center(:,1)=center(:,1)+alignment(f,1); %**ALIGNMENT**
    center(:,2)=center(:,2)+alignment(f,2); %**ALIGNMENT**
    
    [V2cyto, ~, ~]=regionQM(cytolabel, realV2, lo, hi);
    [V2nu, ~, ~]=regionQM(nucmaskV1, realV2, lo, hi);
    [V3nu, ~, ~]=regionQM(nucmaskV1, realV3, lo, hi);

    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numcells=numel(area);
    celldata(1:numcells,f,:)=[center(:,1),center(:,2), area, V2cyto, V2nu, V3nu];
    toc
end

save([datadir,'celldata.mat'],'celldata');
%% track; each individual cell is one layer

maxjump=50;
areadiffthreshold=300; %**AREA SCREEN**
%%% Initialize arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=find(~isnan(celldata(:,1,1)),1,'last');
totalframes=size(celldata,2);
totalfeatures=size(celldata,3);
tracedata=ones(numcells,totalframes,totalfeatures)*NaN;
tracepath=ones(numcells,totalframes)*NaN;

for i=1:numcells
    tracedata(i,1,:)=celldata(i,1,:);
    tracepath(i,1)=i;
    prevID=i;
    xprev=celldata(prevID,1,1);
    yprev=celldata(prevID,1,2);
    areaprev=celldata(prevID,1,3);
    for f=2:totalframes
        %%% update current frame data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xcur=celldata(:,f,1);
        ycur=celldata(:,f,2);
        areacur=celldata(:,f,3);
        %%% find current cells within maxjump from previous frame %%%%%%%%%
        neighbors=find(abs(xcur-xprev)<maxjump & abs(ycur-yprev)<maxjump);
        %%% if neighbors, stop tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(neighbors)
            break;
        end
        %%% find closest cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dist=sqrt((xcur(neighbors)-xprev).^2+(ycur(neighbors)-yprev).^2);
        distidx=find(min(dist));
        match=neighbors(distidx);
        %%% **AREA SCREEN** %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        areadiff=abs(areacur(match)-areaprev);
        if areadiff<areadiffthreshold
            tracedata(i,f,:)=celldata(match,f,:);
            tracepath(i,f)=match;
            xprev=xcur(match);
            yprev=ycur(match);
            areaprev=areacur(match);
            prevID=match;
        else
            break;
        end
    end
end

%% Detect tracking conflicts and remove data after the conflict %%%%%%%%%%

for f=2:totalframes
    % for each frame, find cells redundantly claimed by different tracks.
    sortedcells=sort(tracepath(:,f));
    redundantsortedindex=find(diff(sortedcells)==0);
    % Y = diff(X) calculates differences between adjacent elements of X along the first array 
    conflictcells=unique(sortedcells(redundantsortedindex));
    % for tracks that claim the same cell in a given frame, delete the data
    % for both tracks for that frame and all subsequent frames.
    if ~isempty(conflictcells)
        conflictingtracks=find(ismember(tracepath(:,f),conflictcells));
        tracedata(conflictingtracks,f:totalframes,:)=NaN;
        tracepath(conflictingtracks,f:totalframes)=NaN;
    end
end

save([datadir,'tracedata.mat'],'tracedata');

%% pick only successfully traced cells

tracesuccessindex = ~isnan(tracedata(:, totalframes,1));
tracessid = find(tracesuccessindex>0); %row num for tracepath
tracesuccesscount = length(tracessid);

tracesuccess = tracedata(tracessid',:,:);

save([datadir,'tracesuccess.mat'],'tracesuccess');

%% calculate V2 ratio and rise time, calculate V3 rise time, calculate G1 length
% tracesuccess[center(:,1),center(:,2), area, V2cyto, V2nu, V3nu];

V2ratio = tracesuccess(:, :, 4)./tracesuccess(:, :, 5);
V3nu = tracesuccess(:, :, 6);

V2time = zeros(tracesuccesscount, 1);
V3time = zeros(tracesuccesscount, 1);

for ts = 1:tracesuccesscount
    V2smooth = smooth(V2ratio(ts, :), 20);
    V3smooth = smooth(V3nu(ts, :), 20);
    
    V2diff = smooth(diff(V2smooth), 20);
    V3diff = smooth(diff(V3smooth), 20);
    
    V2T = find(V2diff > quantile(V2diff, 0.6));
    V3diffT = V3diff(V2T:totalframes);
    V3T = find(V3diffT > quantile(V3diffT, 0.6));
    
    V2time(ts) = V2T;
    V3time(ts) = V3T;
end

G1length = V3time - V2time;
