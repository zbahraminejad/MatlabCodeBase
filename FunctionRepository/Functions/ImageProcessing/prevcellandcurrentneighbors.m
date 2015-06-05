function prevcellandcurrentneighborsIF(debugpackage,nuc_mask,winrad,xprev,yprev,previd)
nuc_mask_prev=debugpackage{1};
prevjitter=debugpackage{2}; prevjitx=prevjitter(1); prevjity=prevjitter(2);
reljitter=debugpackage{3}; reljitx=reljitter(1); reljity=reljitter(2);

xprev=xprev-prevjitx; yprev=yprev-prevjity;
fprintf('x=%0.0f y=%0.0f\n',xprev(previd),yprev(previd));

[height,width]=size(nuc_mask);
dxminprev=max([round(xprev(previd)-winrad) 1]); dxmaxprev=min([round(xprev(previd)+winrad) width]);
dyminprev=max([round(yprev(previd)-winrad) 1]); dymaxprev=min([round(yprev(previd)+winrad) height]);
dxmincur=round(dxminprev-reljitx); dxmindiff=double((1-dxmincur)*(dxmincur<1)); dxmincur=max([dxmincur 1]);
dxmaxcur=round(dxmaxprev-reljitx); dxmaxdiff=double((dxmaxcur-width)*(dxmaxcur>width)); dxmaxcur=min([dxmaxcur width]);
dymincur=round(dyminprev-reljity); dymindiff=double((1-dymincur)*(dymincur<1)); dymincur=max([dymincur 1]);
dymaxcur=round(dymaxprev-reljity); dymaxdiff=double((dymaxcur-height)*(dymaxcur>height)); dymaxcur=min([dymaxcur height]);
dbmaskcur=bwmorph(nuc_mask,'remove');
dbmaskcur=dbmaskcur(dymincur:dymaxcur,dxmincur:dxmaxcur);
dbmaskcur=padarray(dbmaskcur,[dymindiff dxmindiff],'pre');
dbmaskcur=padarray(dbmaskcur,[dymaxdiff dxmaxdiff],'post');
dbmaskprev=bwmorph(nuc_mask_prev,'remove');
dbmaskprev=dbmaskprev(dyminprev:dymaxprev,dxminprev:dxmaxprev);
dbimage=mat2gray(dbmaskprev);
dbimage(:,:,2)=dbmaskcur;
dbimage(:,:,3)=0;
figure,imshow(imresize(dbimage,3));