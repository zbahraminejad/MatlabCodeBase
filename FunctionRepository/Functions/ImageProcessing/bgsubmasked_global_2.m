function real=bgsubmasked_global_2(raw,nanmask,numblocks,compression,sampleprctile)
rawcomp=imresize(raw,1/compression,'nearest');
nanmaskcomp=imresize(nanmask,1/compression,'nearest');

fg=rawcomp; fg(nanmaskcomp)=NaN;
tfg=fg; tfg(nanmaskcomp)=[];
vals=tfg(:);
binmax=prctile(vals,95);
binmin=prctile(vals,5);
vals=vals(vals<binmax & vals>binmin);
[kval,xval]=ksdensity(vals);
globalbg=xval(find(kval==max(kval),1,'first'));
if numblocks>1
    [height, width]=size(fg);
%     bg=blocksmooth_mode_3(fg,numblocks,globalbg);
    blockimg=blockpercentile_blockimage(fg,numblocks,sampleprctile);
    blockimg(isnan(blockimg))=globalbg;
    bg=imresize(blockimg,[height width],'bicubic');
    bg=imresize(bg,compression,'bicubic');
else
    bg=ones(size(raw))*globalbg;
end
real=raw-bg;