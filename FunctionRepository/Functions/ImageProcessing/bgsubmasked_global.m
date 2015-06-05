function real=bgsubmasked_global(raw,nanmask,numblocks,compression)
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
    bg=blocksmooth_mode_3(fg,numblocks,globalbg);
    bg=imresize(bg,compression,'bicubic');
else
    bg=ones(size(raw))*globalbg;
end
real=raw-bg;