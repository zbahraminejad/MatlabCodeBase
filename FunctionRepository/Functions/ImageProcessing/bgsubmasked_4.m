function real=bgsubmasked_4(raw,nanmask,numblocks,compression)
rawcomp=imresize(raw,1/compression,'nearest');
nanmaskcomp=imresize(nanmask,1/compression,'nearest');

fg=rawcomp; fg(nanmaskcomp)=NaN;
vals=raw(:);
binmax=prctile(vals,95);
binmin=prctile(vals,5);
vals=vals(vals<binmax & vals>binmin);
[kval,xval]=ksdensity(vals);
globalbg=xval(kval==max(kval));
bg=blocksmooth_mode_3(fg,numblocks,globalbg);
bg=imresize(bg,compression,'bicubic');
real=raw-bg;