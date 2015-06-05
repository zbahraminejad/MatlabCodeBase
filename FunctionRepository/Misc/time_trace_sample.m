%%% Function To Make Sequence of Frames from input time trace
row=3;col=2;site=1;
cellID = 2528;

imagepath = 'G:\Michael\';
experimentpath='20150401-CC-Diff\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
wellname = nameandsite(shot);
rawdir=[imagepath,experimentpath,'Raw\',wellname,shot,'_'];
datadir=([imagepath,experimentpath,'Data\']);
datafile = [datadir,'tracedata_',shot,'_nolink','.mat'];
motheroption = 0;
daughteroption =0;
IFoption = 0;
load(datafile,'tracedata','genealogy','jitters');
[tracedata,tracestats,motherstats,IFdata,IDs,markedmitosis,lastcellmother]=gathertracedata_mz_1(datadir,datafile,shot,motheroption,daughteroption,IFoption);
figure,
plot(tracedata(cellID,:,10));
frames  = 1:24:653;
numframes = length(frames);
nucr = 8;
windowsize = 2*nucr;

channelNames = {'Cy5_','YFP_','mCherry_'};
numchannels = numel(channelNames);
testimage = imread([rawdir,channelNames{1},num2str(frames(1)),'.tif']);
[imheight,imwidth] = size(testimage);
 
figure('units','normalized','outerposition',[0 0 1 1]);
for fcount = 1:numframes
    xcent = int16(tracedata(cellID,frames(fcount),1)-jitters(frames(fcount),1));
    ycent = int16(tracedata(cellID,frames(fcount),2)-jitters(frames(fcount),2));
    xmin = max([xcent-windowsize 1]);
    xmax = min([xcent+windowsize imwidth]);
    ymin = max([ycent-windowsize 1]);
    ymax = min([ycent+windowsize imheight]);
    
    currimage1 = double(imread([rawdir,channelNames{1},num2str(frames(fcount)),'.tif']));
    currimage2 = double(imread([rawdir,channelNames{2},num2str(frames(fcount)),'.tif']));
    currimage3 = double(imread([rawdir,channelNames{3},num2str(frames(fcount)),'.tif']));
    croppedcurr1 = currimage1(ymin:ymax,xmin:xmax);
    croppedcurr2 = currimage2(ymin:ymax,xmin:xmax);
    croppedcurr3 = currimage3(ymin:ymax,xmin:xmax);

    
    R = zeros(length(ymin:ymax),length(xmin:xmax),3); 
    G = zeros(length(ymin:ymax),length(xmin:xmax),3); 
    B = zeros(length(ymin:ymax),length(xmin:xmax),3); 
    
    R(:,:,1) = mat2gray(croppedcurr3);
    G(:,:,2) = mat2gray(croppedcurr2);
    B(:,:,3) = mat2gray(croppedcurr1);
    
    subtightplot(numchannels,numframes,fcount,[0 0],0.4,0.05), imshow(B,[0 1000]);
    subtightplot(numchannels,numframes,fcount+numframes,[0 0],0.4,0.05),imshow(G,[0 1000]);
    subtightplot(numchannels,numframes,fcount+numframes*2,[0 0],0.4,.05),imshow(R,[0 1000]);
    
    
end



%%%%%Mark the spot along the trace where we sampled from%%%%%%%%%%%%%%%%%
% timetrace1 = traces1(1,:);
% timetrace2 = traces2(1,:);

