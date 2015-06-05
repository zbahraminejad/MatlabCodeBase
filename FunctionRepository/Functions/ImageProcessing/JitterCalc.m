function [xcoor,ycoor,TempStatM,minV]=JitterCalc(ImageO,ImageN,x0,y0,r)

%ImageO:original,ImageN:new,[x0,y0]:starting coor,r:grid
%[xcoor,ycoor]=resulting coor
%ImageN=YfpImage;ImageO=YfpImage;r=4;x0=4;y0=4;

%% preparation
ymatrix=y0+(-3:3)'*ones(1,7)*r;xmatrix=x0+ones(7,1)*(-3:3)*r;

%% calculation
TempStat=zeros(49,3);
TempStatM=zeros(7,7);
for rows=1:7
    for cols=1:7

        dy=ymatrix(rows,cols);dx=xmatrix(rows,cols);

        ImageN2=ImageN(max(1,1-dy):min(end,end-dy),max(1,1-dx):min(end,end-dx));
        ImageO2=ImageO(max(1,1+dy):min(end,end+dy),max(1,1+dx):min(end,end+dx));
        ImageD=abs(ImageN2-ImageO2);
        
        TempStat(rows*7-7+cols,:)=[dx,dy,sum(ImageD(:))/size(ImageD(:),1)];
        TempStatM(rows,cols)=sum(ImageD(:))/size(ImageD(:),1);

    end
end

%% sorting and get resulting coor
[minV,I]=min(TempStat(:,3));
xcoor=TempStat(I,1);ycoor=TempStat(I,2);