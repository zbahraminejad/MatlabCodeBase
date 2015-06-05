%function moviemaker_test(row,col,site)
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
label = 1;
%nucroot = '_nucedge_';
nucroot = '_nucedge&ring_';
hDHBroot = '_hDHB_';
CEroot = '_cdt1_';

%%% Prepare Movie
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
M=VideoWriter([datadir,shot,'_numvid_outlier_march'],'Uncompressed AVI');
M.FrameRate = 4;
open(M);

%%% Retrieve data for cell ID and locations
if label
    load([datadir,shot,'_alldata_outlier_march'],'bestsp','best_rc','corruptlist','leaveoutlist')
    totalcells = size(bestsp{end},1);
    totalframes = size(bestsp,3);
    cellidlist = [1:totalcells];
    leavein = cellidlist(~ismember(cellidlist,[leaveoutlist;corruptlist]));
end
%totalframes = 5;
ra=1;
%for f=1:totalframes
for f=1
    rawimage = single(imread([cpdir,shot,nucroot,num2str(f),'.tif']));
    rawimage2 = single(imread([cpdir,shot,hDHBroot,num2str(f),'.tif']));
    tempframe = imadjust(mat2gray(rawimage));
    tempframe(:,:,2) = imadjust(mat2gray(rawimage2));
    tempframe(:,:,3) = zeros(size(rawimage));
    
    if label
        possiblecells=find((bestsp{f}(:,1)>0));
        actualcells=possiblecells(ismember(possiblecells,leavein));
        imshow(tempframe)
        title(['Frame ',num2str(f)],'fontweight','bold')
        hold on;
        for cc=1:length(possiblecells)
            i=possiblecells(cc);
            text(ra*(bestsp{f}(i,1)),ra*(bestsp{f}(i,2))+round(6*ra),num2str(i),'horizontalalignment','center','color','w','fontsize',12,'fontweight','bold')
        end
        writeVideo(M,getframe(gcf));
    else
        writeVideo(M,im2frame(tempframe));
    end
end
close(M);