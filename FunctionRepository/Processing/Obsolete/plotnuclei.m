tempdir='h:\Documents\Timescape\20120807_12drugs\';  %change this
subdir='Data\';

for row=3
    for col=2
        for site=0  

tempsite=[num2str(row),'_', num2str(col), '_', num2str(site)];


load([tempdir,subdir,tempsite,'_alldata.mat'])
movienonums=aviread([tempdir,subdir,tempsite,'.avi']);


ce='rgbmkw';cf='grycwk';
ra=1;

%%
moviename= strcat(tempsite,'_numberedMovie.avi');
M=avifile([tempdir,subdir,moviename],'compression','none','fps',4);
%M=avifile([tempdir,subdir,'tempsite','/MovieWithNuclei_New.avi'],'compression','none','fps',4);

figure;hold off;set(gcf,'color','w')
for f=1:size(bestsp,3)
    disp(f)
    %%

    tempframe=movienonums(f).cdata;   %contains the 3 layers of color
    
    imshow(tempframe)
    title(['Frame ',num2str(f)],'fontweight','bold')
    
    %% generate nuclei info
    lo=find((bestsp{f}(:,1)>0));
    hold on
    for cc=1:length(lo)
        cn=mod(lo(cc),6)+1;
        plot(ra*(bestsp{f}(lo(cc),1)),ra*(bestsp{f}(lo(cc),2)),'o','markersize',round(3*ra),'markeredgecolor',ce(cn),'markerfacecolor',cf(cn))
        text(ra*(bestsp{f}(lo(cc),1)),ra*(bestsp{f}(lo(cc),2))+round(6*ra),num2str(lo(cc)),'horizontalalignment','center','color','k','fontsize',10,'fontweight','bold')
    end
    hold off
    %%
    M=addframe(M,getframe(gcf));
end

%%
M=close(M);
close(gcf)

        end
    end
end