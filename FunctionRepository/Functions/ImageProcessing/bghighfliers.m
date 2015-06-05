function [real,highfliers]=bghighfliers(raw,nanmask,nucr,iqrmult,numblocks)
blur=imfilter(raw,fspecial('disk',nucr),'symmetric');
highfliers=maskhighfliers_nonnuc(blur,nanmask,nucr,iqrmult);
fg=nanmask | highfliers;
blur(fg)=NaN;
bg=blocksmooth_mode(blur,numblocks);
real=raw-bg;