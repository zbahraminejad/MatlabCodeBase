function cmap=makecmap(firstcode,midcode,lastcode,midpoint)
cmap = ones(128,3)*NaN;
firststep=(midcode-firstcode)/(midpoint-1);
for i=1:3
    if firststep(i)==0
        panel=ones(1,midpoint)*firstcode(i);
    else
        panel=firstcode(i):firststep(i):midcode(i);
    end
    cmap(1:midpoint,i)=panel;
end
secondpointnum=128-midpoint+1;
secondstep=(lastcode-midcode)/(secondpointnum-1);
for i=1:3
    if secondstep(i)==0
        panel=ones(1,secondpointnum)*midcode(i);
    else
        panel=midcode(i):secondstep(i):lastcode(i);
    end
    cmap(midpoint:128,i)=panel;
end