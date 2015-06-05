function numcycling = countcycling(traces,tracestats)
[numcells,~] = size(traces);
tempindex = zeros(numcells,1);
for j = 1:numcells
   findchange = traces(j,tracestats(j,1):end)> 1;
   if sum(findchange)<12
       tempindex(j) = -1;
   else
       tempindex(j) = find(findchange,1,'first');
   end
   
end
numcycling = sum(tempindex>0);