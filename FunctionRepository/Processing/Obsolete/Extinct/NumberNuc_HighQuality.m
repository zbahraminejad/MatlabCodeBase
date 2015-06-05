%%
%Mingyu
clear mex;
path = 'h:\Documents\Timescape\20120807_12drugs\';  %folder containing the movies folders [CHANGE]
datadir = [path,'DataEval\'];
rawdir = [path,'Raw\'];
cpdir = [path,'CroppedProcessed\'];
SF=100;
EF=105;

for row=3
    for col=1
        for site=1  

tempsite=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,tempsite,'_alldata.mat'])

ce='rgbmkw';cf='grycwk';
ra=1;

%%
moviename= strcat(tempsite,'_numberedMovie.avi');
M=avifile([datadir,moviename],'compression','none','fps',4);

figure;hold off;%set(gcf,'color','w')
for f=1:size(bestsp,3)
    disp(f)
    %%
    frame=f+SF-1;
    DAs_or=single(imread([cpdir,movieName,'_mRFP1_',num2str(frame),'.tif']));
    REs_or=single(imread([cpdir,movieName,'_EYFP_',num2str(frame),'.tif']));
    CEs_or=single(imread([cpdir,movieName,'_ECFP_',num2str(frame),'.tif']));
    %% Add frame to movie
    tempframe=imadjust(mat2gray(DAs_or)); %red
    tempframe(:,:,2)=imadjust(mat2gray(REs_or)); %green.  
    tempframe(:,:,3)=imadjust(mat2gray(CEs_or)); %blue
    
    imshow(tempframe);
    title(['Frame ',num2str(frame)],'fontweight','bold');
    
    %% generate nuclei info
    lo=find((bestsp{f}(:,1)>0));
    hold on
    for cc=1:length(lo)
        cn=mod(lo(cc),6)+1;
        %plot(ra*(bestsp{f}(lo(cc),1)),ra*(bestsp{f}(lo(cc),2)),'o','markersize',round(3*ra),'markeredgecolor',ce(cn),'markerfacecolor',cf(cn))
        text(ra*(bestsp{f}(lo(cc),1)),ra*(bestsp{f}(lo(cc),2))+round(6*ra),num2str(lo(cc)),'horizontalalignment','center','color','k','fontsize',10,'fontweight','bold');
    end
    hold off
    %%
    M=addframe(M,im2frame(tempframe));
end

%%
M=close(M);
%close(gcf)

        end
    end
end