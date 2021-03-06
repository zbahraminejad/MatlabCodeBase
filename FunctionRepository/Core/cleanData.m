function [meansite,datadist] = cleanData(data,site_num)
    % data = percentile5to95(data);

	%---------------------------------------------
	% i = 1: data(1:site_num)
	% i = 2: data(site_num + 1:2 * site_num)
	% ...
	% [1:site_num, site_num + 1:2 * site_num, ...]
	%---------------------------------------------
    for i = 1:length(data)/site_num
		% Convert cell array to ordinary array of the underlying data type
        % Partition the data array based on the specified number of sites
        datadist{i} = cell2mat(data((i - 1) * site_num + 1:site_num * i));
    end

    % The input data array comprises a vertical array of doubles
    % Apply 'mean' function upon each cell
    meansite = cellfun(@mean,data);

    if site_num ~= 1
		% Reshape the array of data into a matrix
        meansite = reshape(meansite,site_num,length(meansite)/site_num);
        meansite = mean(meansite);
    end
end

function [store] = percentile5to95(data)
	store = cell(length(data),1);
	for i = 1:length(data)
    	bottom = prctile(data{i},5);
    	top = prctile(data{i},95);
    	store{i} = data{i}(data{i} < top & data{i} > bottom);
	end
end