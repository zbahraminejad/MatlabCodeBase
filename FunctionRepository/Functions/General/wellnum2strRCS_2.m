function [rowstring,colstring,sitestring] = wellnum2str_2(row,col,site)
rowstringtable = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
colstringtable = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
sitestringtable = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36'};
rowstring = rowstringtable{row};
colstring = colstringtable{col};
sitestring = sitestringtable{site};
end