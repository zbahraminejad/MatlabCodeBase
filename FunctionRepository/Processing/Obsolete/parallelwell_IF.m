format short;
cd Functions
rows = [2 3 4 5 7];
numrows=length(rows);
parfor idx=1:numrows
    row=rows(idx);
    fprintf('Row %0.0f\n',row);
    X_immunofluorescence_simple(row);
end