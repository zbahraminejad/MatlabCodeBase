function real=bgsubmasked(raw,nanmask,numblocks)
fg=raw; fg(nanmask)=NaN;
%bg=blocksmooth(fg,10);
bg=blocksmooth_mode_2(fg,numblocks);
real=raw-bg;