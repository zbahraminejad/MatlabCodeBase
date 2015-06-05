function mtx=ppdist2(mt1,mt2)

%%
dm1=size(mt1,2);dm2=size(mt2,2);
sz1=size(mt1,1);sz2=size(mt2,1);
pt1=ceil(sz1/1000);pt2=ceil(sz2/1000);

%%
if dm1==dm2
    mtx=zeros(sz1,sz2);
    for row=1:pt1
        for col=1:pt2
            x1=1+1000*(col-1);x2=min([1000*col,sz2]);
            y1=1+1000*(row-1);y2=min([1000*row,sz1]);
            for dd=1:dm1
                if dd==1
                    tmt=zeros(y2-y1+1,x2-x1+1);
                end
                tmt=tmt+(ones(y2-y1+1,1)*(mt2(x1:x2,dd)')-mt1(y1:y2,dd)*ones(1,x2-x1+1)).^2;
            end
            mtx(y1:y2,x1:x2)=sqrt(tmt);
        end
    end
end