function real=bgsubmasked_3(raw,nanmask,numblocks)
fg=raw; fg(nanmask)=NaN;
vals=raw(:);
binmax=prctile(vals,95);
binmin=prctile(vals,5);
vals=vals(vals<binmax & vals>binmin);
[kval,xval]=ksdensity(vals);
globalbg=xval(kval==max(kval));
bg=blocksmooth_mode_3(fg,numblocks,globalbg);
real=raw-bg;