function colorcode=distributecolormap(colorkey,condnum)
cmap=colormap(colorkey); close(gcf);
colorcode=ones(condnum,3);
colorcode(1,:)=cmap(1,:); colorcode(condnum,:)=cmap(64,:);
cmapstep=64/(condnum-1);
for i=1:condnum-2
    colorcode(i+1,:)=cmap(round(cmapstep*i),:);
end