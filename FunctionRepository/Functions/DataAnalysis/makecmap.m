function cmap=makecmap(lowcode,midcode,highcode,anchorlow,anchormid,anchorhigh)
cmap = ones(64,3)*NaN;
firststep=(midcode-lowcode)/(anchormid-anchorlow); %will set the first 32 or 33 values
for i=1:3 %each color separately
    if firststep(i)==0
        panel=ones(1,anchormid)*lowcode(i);
    else
        panelanchorlowtomid=lowcode(i):firststep(i):midcode(i);
        if anchorlow>1
            panelmintoanchorlow=ones(1,anchorlow-1)*lowcode(i);
            panel=[panelmintoanchorlow,panelanchorlowtomid];
        else
            panel=panelanchorlowtomid;
        end
    end
    cmap(1:anchormid,i)=panel;
end
secondpointnum=anchorhigh-anchormid+1; %will set the last 33 values or 32 values (will overwrite midval w/ same val)
secondstep=(highcode-midcode)/(secondpointnum-1);
secondtotalpointnum=64-anchormid+1;
for i=1:3
    if secondstep(i)==0
        panel=ones(1,secondtotalpointnum)*midcode(i);
    else
        panelanchormidtohigh=midcode(i):secondstep(i):highcode(i);
        if anchorhigh<64
            panelanchorhightomax=ones(1,64-anchorhigh)*highcode(i);
            panel=[panelanchormidtohigh,panelanchorhightomax];
        else
            panel=panelanchormidtohigh;
        end
    end
    cmap(anchormid:64,i)=panel;
end