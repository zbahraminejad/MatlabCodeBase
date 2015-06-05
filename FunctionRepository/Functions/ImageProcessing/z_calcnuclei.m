%% preparation
clc;close all;clear all
time1=tic;
format short
ry1=20;ry2=20;rx1=50;rx2=50;
tempdir='F:\AxonDataBase\20110310_drug';
mkdir([tempdir,'/data'])
SF=1;EF=49;%SF=starting frame, EF = ending frame
%sss=cell(8,12,EF); % used for storing cell coordinates
nucr=3;

%% loop
for i=2:7 %row indecies
    alp=char(64+i);
    for j=3:10 %column indecies
        %% setup tempwell
        if j<10
            tempwell=[alp,'0',num2str(j)];%alp= alphabetic letter, '0' , num2str fxn assigns a character to the value
        else tempwell=[alp,num2str(j)];
        end
        x=0;y=0;%store jitters data
        welldir=dir([tempdir,'/',tempwell,'*.tif']);
        wellsss=cell(1,1,EF-SF+1);
        %%
        M=avifile([tempdir,'/data/',tempwell,'.avi'],'fps',4,'compression','Cinepak');
        
        %% dapi
        for f=SF:EF
            disp([i,j,f])%show row, column, frame
            time2=tic;
            
            %% reading images
            DAs_or=single(imread([tempdir,'/',welldir(f).name]));%DA = dapi, reads image
            
            %% correcting jitters
            if f>SF
                DI1=imadjust(mat2gray(log(DAs_or)));
                [xx,yy]=CalcJitter(DI0,DI1);
                %if abs(xx)<0.5;xx=0;end
                %if abs(yy)<0.5;yy=0;end
                x=x+xx;y=y+yy;
            end
            DI0=imadjust(mat2gray(log(DAs_or)));
            disp([x,y])
            x0=round(x);y0=round(y);
            ts=size(DAs_or);%size of image
            %x0=0;y0=0;
            
            %% cropping
            DAs_or=DAs_or(1+ry1-y0:ts(1)-ry2-y0,1+rx1-x0:ts(2)-rx2-x0);
            
            %% image processsing
            DAs_bl=imfilter(log(DAs_or),fspecial('disk',round(nucr/1)),'symmetric');
            DAs_bs=bgsub(DAs_bl,3*nucr,0);
            
            %% add frame
            tempframe=imadjust(mat2gray(imresize(DAs_bs,0.5)));
            tempframe(:,:,2)=tempframe;
            tempframe(:,:,3)=tempframe(:,:,1);
            M=addframe(M,im2frame(tempframe));
            
            %% get data
            DAs_ma=getdapimask(DAs_bs,nucr);
            %DAs_la=bwlabel(DAs_ma);
            DAs_da=regionprops(DAs_ma,'Centroid','PixelList','PixelIdxList','Area','Perimeter','EquivDiameter');
            XX=zeros(size(DAs_da,1),1);YY=zeros(size(DAs_da,1),1);
            AC=zeros(size(DAs_da,1),1);%PP=zeros(size(DAs_da,1),1);DI=zeros(size(DAs_da,1),1);
            DD=zeros(size(DAs_da,1),1);
            for cc=1:size(DAs_da,1)
                XX(cc,1)=DAs_da(cc).Centroid(1);
                YY(cc,1)=DAs_da(cc).Centroid(2);
                AC(cc,1)=DAs_da(cc).Area;
                %                 PP(cc,1)=DAs_da(cc).Perimeter;
                %                 DI(cc,1)=DAs_da(cc).EquivDiameter;
                DD(cc,1)=mean(DAs_bs(DAs_da(cc).PixelIdxList));
            end
            
            %% filter data
%             HT=max([400,10*mean(AC)]);%mean(AC0(AC0>LT))+3.29*std(AC0(AC0>LT));
%             if f==SF
%                 [~,LT0]=ThreshImage_smbg(AC);
%                 [~,~,DLT0]=ThreshImage_smbg(DD);
%                 range0a=median(AC(AC>LT0));
%                 %range0d=diff(prctile(DD,[0.1,99.9]));
%                 range0d=median(DD(DD>DLT0));
%             end
%             %disp([LT,HT])
%             range1d=median(DD(DD>DLT0));
%             %range1d=diff(prctile(DD,[0.1,99.9]));
%             range1a=median(AC(AC>LT0));
%             LT=LT0*(range1a/range0a);
%             DLT=DLT0*(range1d/range0d);
%             lo=((DD>=DLT)&(AC>=LT))&(AC<=HT);
%             wellsss{:,:,f-SF+1}=[XX(lo),YY(lo),DD(lo),AC(lo)];
            wellsss{:,:,f-SF+1}=[XX,YY,DD,AC];
            toc(time2)
        end
        save([tempdir,'/data/well_',num2str(i),'_',num2str(j),'.mat'],'wellsss')
        %%
        M=close(M);
    end
end

%% end
toc(time1)