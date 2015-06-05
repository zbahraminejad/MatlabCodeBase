function [xx,yy,minV]=CalcJitter(DapiImage0,DapiImage)
% DapiImage0:original image, DapiImage:new image

%% Get mask
%Mask0=GetM(DapiImage0);Mask=GetM(DapiImage);
Mask0=DapiImage0;Mask=DapiImage;

%% Calculation
xcoor=0;ycoor=0;
[sy,sx]=size(DapiImage);

for r=[64,16,4,1]   %r=[16,4,1; r=[4,1]  %first sweep with stepsize of 4 pixels, then sweep with step size of 1 pixel. [4,1] is the original.
    x0=xcoor+1;y0=ycoor+1;%make while loop running
    while ((xcoor~=x0)||(ycoor~=y0))&&((abs(xcoor)<sx)&&(abs(ycoor)<sy))
        x0=xcoor;y0=ycoor;
        [xcoor,ycoor,TempStatM,minV]=JitterCalc(Mask0,Mask,x0,y0,r);
    end
end

%% get better value

px=zeros(1,3);py=zeros(1,3);
ts=TempStatM;
for rd=1:3
    ts2=imresize(ts,[61,61],'bicubic');
    ts3=ts2;  ts3(:,1:3)=[]; ts3(:,end-2:end)=[]; ts3(1:3,:)=[]; ts3(end-2:end, :)=[]; %take out the outer 3 rows and 3 columns    
    [dummy,tpx]=min(min(ts3,[],1));
    [dummy,tpy]=min(min(ts3,[],2));
    tpx=tpx+3; tpy=tpy+3;
    px(1,rd)=(tpx-31)/(10^rd);
    py(1,rd)=(tpy-31)/(10^rd);   
    ts=ts2((tpy-3):(tpy+3),(tpx-3):(tpx+3));
end
xx=xcoor+sum(px);
yy=ycoor+sum(py);

% TempStatM2=TempStatM(3:5,3:5);
% TempStatM3=imresize(TempStatM2,[201,201],'bicubic');
% [~,px]=min(min(TempStatM3,[],1));
% [~,py]=min(min(TempStatM3,[],2));
% TempStatM4=imresize(TempStatM3(py-1:py+1,px-1:px+1),[201,201],'bicubic');
% [~,px2]=min(min(TempStatM4,[],1));
% [~,py2]=min(min(TempStatM4,[],2));
% xx=xcoor+(px-101)/100+(px2-101)/10000;
% yy=ycoor+(py-101)/100+(py2-101)/10000;
