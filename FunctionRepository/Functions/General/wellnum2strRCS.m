function [rowstring,colstring,sitestring] = wellnum2str(row,col,site)
rowstringtable = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
colstringtable = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
sitestringtable = {'1' '2' '3' '4' '5' '6' '7' '8' '9'};
rowstring = rowstringtable{row};
colstring = colstringtable{col};
sitestring = sitestringtable{site};
end