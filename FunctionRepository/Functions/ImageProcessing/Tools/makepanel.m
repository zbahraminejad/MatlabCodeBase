function makepanel(nuc_raw,mask,centers,winrad,panelsperside)
%gate out centers (x,y of each cell) to population of interest ahead of time.
%winrad usually 3*nucr.  panelsperside usually 5
x=centers(:,1); y=centers(:,2);
[height width]=size(nuc_raw);
winlength=2*winrad+1;
inbounds=find(x>winrad+1 & x<width-winrad & y>winrad+1 & y<height-winrad);
x=x(inbounds); y=y(inbounds);

minx=ceil(x-winrad); maxx=minx+winlength-1;
miny=ceil(y-winrad); maxy=miny+winlength-1;
%%% make panel border template %%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempoutline=zeros(winlength,winlength);
tempoutline(1:end,[1 end])=1; tempoutline([1 end],1:end)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numpanels=numel(inbounds);
minint=ones(numpanels,1)*NaN;
for k=1:numpanels
    j=mod(k,panelsperside^2); j=j+panelsperside^2*(j==0);
    panelx=mod(k,panelsperside); panelx=panelx+panelsperside*(panelx==0); panelx=panelx-1;
    panely=ceil(j/panelsperside)-1;
    panelx=panelx*winlength+1;
    panely=panely*winlength+1;
    windowraw=nuc_raw(miny(k):maxy(k),minx(k):maxx(k));
    minint(j)=prctile(windowraw(:),10);
    panelraw(panely:panely+winlength-1,panelx:panelx+winlength-1)=windowraw;
    panelmask(panely:panely+winlength-1,panelx:panelx+winlength-1)=mask(miny(k):maxy(k),minx(k):maxx(k));
    paneloutline(panely:panely+winlength-1,panelx:panelx+winlength-1)=tempoutline;
    if j==panelsperside^2 || k==numpanels
        panelraw(panelraw==0)=nanmean(minint);
        tempframe=imadjust(mat2gray(panelraw));
        tempframe(:,:,2)=panelmask;
        tempframe(:,:,3)=paneloutline;
        figure,imshow(tempframe);
        panelraw=0; panelmask=0; paneloutline=0;
    end
end