function real=subtractbg(raw,nanmask)
rawnan=raw;
rawnan(nanmask)=NaN;
globalbg=nanmedian(rawnan(:));
real=raw-globalbg;
end