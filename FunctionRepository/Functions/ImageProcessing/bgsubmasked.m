function real=bgsubmasked(raw,nanmask,numblocks)
fg=raw; fg(nanmask)=NaN;
%bg=blocksmooth(fg,10);
bg=blocksmooth_mode(fg,numblocks);
real=raw-bg;