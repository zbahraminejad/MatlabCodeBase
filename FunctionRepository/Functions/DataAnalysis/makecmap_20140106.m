function cmap=makecmap(firstcode,midcode,lastcode)
cmap = ones(128,3)*NaN;
step=(midcode-firstcode)/63;
for i=1:3
    if step(i)==0
        panel=ones(1,64)*firstcode(i);
    else
        panel=firstcode(i):step(i):midcode(i);
    end
    cmap(1:64,i)=panel;
end
step=(lastcode-midcode)/64;
for i=1:3
    if step(i)==0
        panel=ones(1,65)*midcode(i);
    else
        panel=midcode(i):step(i):lastcode(i);
    end
    cmap(64:128,i)=panel;
end