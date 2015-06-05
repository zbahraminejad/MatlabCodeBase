% Classify points

function [classifiedcells,finaldistances,finalindices,highnum,lownum,weirdnum] = classifycell(centroids,centroidlinearind,highind,lowind,weirdind,cellstoconsider)

cents = centroids;
centroidindex = centroidlinearind;
signum = cellstoconsider;

totalcells = length(cents);
distances = pdist(cents);
distances = squareform(distances);

%memory
finaldists = zeros(signum,totalcells);
finalinds = zeros(signum,totalcells);
% meandists = zeros(totalcells,1);

for cellind = 1:totalcells
    [vals,ind]=sort(distances(:,cellind),'ascend');
    ind(1) = [];vals(1)=[];
    finaldists(:,cellind) = vals(1:signum);
    finalinds(:,cellind) = ind(1:signum);
%     meandists(cellind) = mean(finaldists(:,cellind));    
end

[~,highnum,~] = intersect(centroidindex,highind);
[~,lownum,~] = intersect(centroidindex,lowind);
[~,weirdnum,~] = intersect(centroidindex,weirdind);

classify = zeros(signum,totalcells);
for i = 1:totalcells
    for j = 1:signum
        if any(finalinds(j,i)==highnum)
        classify(j,i) = 1;
        elseif any(finalinds(j,i)==lownum)
            classify(j,i) = 2;
        elseif any(finalinds(j,i)==weirdnum)
            classify(j,i) = 3;
        else 
            classify(j,i) = nan;
        end
    end
end
classifiedcells = classify;
finaldistances = finaldists;
finalindices = finalinds;


end