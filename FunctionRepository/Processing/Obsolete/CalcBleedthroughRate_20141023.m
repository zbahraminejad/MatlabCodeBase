projectpath='D:\Documents\Projects\';
imagepath='D:\Images\';
%imagepath='E:\';
shadingpath='D:\Images\ShadingImages\20140424 DCYTC 10x\';
%shadingpath='D:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='20130607 p21dCy2\20140208 H2B DHB p21dCy1dK\';
datadir=([projectpath,experimentpath,'Data\']);
rawdir=[imagepath,experimentpath,'Raw\'];
biasdir=[imagepath,experimentpath,'Raw\Bias\'];
nucname='H2B';
primaryname='p21dCy1dK';
secondaryname='p21';
nucr=12; debrisarea=200;
separatedirectories=1;

%%% wells %%%%%%%%%%%%%%%%%%%%%%%%
rowmat=[3:6];
colmat=[6];
sitemat=1:4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat);
stain='stain';
totalwells=numrows*numcols*numsites;
bleedthroughrates=zeros(totalwells,2);

load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
indexcount=0;
for r=1:numrows
    row=rowmat(r);
    for c=1:numcols
        col=colmat(c);
        for s=1:numsites
            indexcount=indexcount+1;
            site=sitemat(s);
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            if separatedirectories
                rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
            else
                rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
            end
            load([biasdir,nucname,'_',num2str(site),'.mat']); bias1=bias;
            load([biasdir,primaryname,'_',num2str(site),'.mat']); bias2=bias;
            load([biasdir,secondaryname,'_',num2str(site),'.mat']); bias3=bias;
            %%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nuc_raw=single(imread([rawdir,nucname,'_stain.tif'])); nuc_raw=(nuc_raw-bgcmos)./bias1;
            primary_raw=single(imread([rawdir,primaryname,'_stain.tif'])); primary_raw=(primary_raw-bgcmos)./bias2;
            secondary_raw=single(imread([rawdir,secondaryname,'_stain.tif'])); secondary_raw=(secondary_raw-bgcmos)./bias3;

            %%% marker-based watershed: regionalmax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            blurradius=3;
            nuc_mask=threshmask(raw1,blurradius);
            nuc_mask=markershed(nuc_mask,nucr*2/3);
            foreground=nuc_mask;
            nuc_mask=bwareaopen(nuc_mask,debrisarea);
            nuc_mask=secondthresh(raw1,blurradius,nuc_mask,boulderarea);
            nuc_mask=bwareaopen(nuc_mask,debrisarea);
            nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
            nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
            
            compression=4;
            nanmask=imdilate(nuc_mask,strel('disk',nucr/2));
            primary_blur=imfilter(primary_raw,fspecial('gaussian',5),'symmetric');
            primary_real=bgsubmasked_global(primary_blur,nanmask,1,compression);
            secondary_blur=imfilter(secondary_raw,fspecial('gaussian',5),'symmetric');
            secondary_real=bgsubmasked_global(secondary_blur,nanmask,1,compression);
            
            [~,numcells]=bwlabel(nuc_mask);
            nuc_info=regionprops(nuc_mask,'PixelIdxList');
            nanvec=ones(numcells,1)*NaN;
            primary_sig=nanvec;
            secondary_sig=nanvec;
            for cc=1:numcells
                primary_sig(cc)=mean(primary_real(nuc_info(cc).PixelIdxList));
                secondary_sig(cc)=mean(secondary_real(nuc_info(cc).PixelIdxList));
            end
            p=robustfit(primary_sig,secondary_sig);

            %dscatter(primary_sig,secondary_sig);
            %xmin=min(primary_sig); xmax=max(primary_sig); xstep=(xmax-xmin)/100; xvals=xmin:xstep:xmax;
            %hold on; plot(xvals,p(1)+p(2)*xvals,'r');
            %axis([0 2500 0 250]); xlabel('TxRed RFU'); ylabel('Cy5 RFU');
            %set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig.jpg');
            
            bleedthroughrates(indexcount,:)=[p(1) p(2)]; %store [intercept slope]
        end
    end
end

bleedthroughrate=mean(bleedthroughrates);
save([datadir,'bleedthroughrate_p21dCy1dKtop21.mat'],'bleedthroughrate');

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}