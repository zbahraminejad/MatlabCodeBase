function bg=getbg(im,bm,szo)

%% image to start
bm0=im.*bm;

%% set up mask
szi=round(szo*sqrt(0.5));szd=szo-szi;
cir_o=strel('disk',szo,0);cir_o=getnhood(cir_o);
cir_i=strel('disk',szi,0);cir_i=getnhood(cir_i);
cir_r=cir_o;
cir_r(1+szd:end-szd,1+szd:end-szd)=~cir_i;
cir_r(1+szo,:)=[];cir_r(:,1+szo)=[];
cir_r=cir_r/sum(cir_r(:));

%%
bg=bm0;bm2=bm;
bmnan=1;
while bmnan>0
    bm1=imfilter(bg,cir_r,'symmetric');
    %figure;imagesc(bm1)
    bm2=imfilter(single(bm2>0),cir_r,'symmetric');
    bm3=bm1./bm2;
    bm3nan=isnan(bm3);bm3(bm3nan)=0;
    bmnan=sum(bm3nan(:));
    bg=bm3;
end