function drawHistTaka(conditions,distData,bmin,bmax,varargin)

%---Arguments Arragement----
% first:conditions
% second:distData
% third:bmin
% fourth:bmax
% fifth:numberOfHistogramBins
%---------------------------

    numberOfConditions = size(conditions,1);
    numberOfHistogramBins = 30;
    if size(varargin,2) > 0
        if ~isempty(varargin{1})
            numberOfHistogramBins = varargin{1};
        end
    end

    bstep = (bmax - bmin)/numberOfHistogramBins; % Prob. Dist. Fun. 30
    histogramBin = bmin:bstep:bmax;
    binfill = [histogramBin fliplr(histogramBin)];
    namecount = cell(numberOfConditions,1);
    colorcode = distributecolormap(jet,numberOfConditions); % alt: cool
    colors = colorcode+.01*(1-colorcode);

    figure; hold on;

    for i = 1:numberOfConditions
        rowmat = cell2mat(conditions(i,2));
        colmat = cell2mat(conditions(i,3));
        dat = [];
        for row = rowmat
            for col = colmat
                    dat = [dat; distData{row,col}];
            end
        end

        % bottom = prctile(dat,5);
        % top = prctile(dat,95);
        % dat = dat(dat<top&dat>bottom);

        if size(varargin,2) == 2
            if strcmp(varargin{2},'log') ~= 0
                dat = log(dat);
                dat = abs(dat);
            end
        end

        histData = histc(dat,histogramBin);
        histData = 100*histData/sum(histData);hold on
        fill(binfill,[histData;zeros(length(histData),1)],colors(i,:),...
            'edgecolor',colors(i,:),'linewidth',1.5,'FaceAlpha', 0.4);

        namecount{i} = char(conditions(i,1));
    end

    legend(char(namecount(:)),'location','northeast');

    % xlabel('log2(meanRFU)');
    ylabel('pdf (%)'); 
    xlim([bmin bmax]);
    set(gcf,'color','w','PaperPosition',[0 0 8 6]);
end

function [store] = percentile5to95(data)

    store = cell(length(data),1);
    for i = 1:length(data)
        bottom = prctile(data{i},5);
        top = prctile(data{i},95);
        store{i} = data{i}(data{i} < top & data{i} > bottom);
    end
end

function colorcode = distributecolormap(colorkey,numberOfConditions)
    cmap = colormap(colorkey); close(gcf);
    colorcode = ones(numberOfConditions,3);
    colorcode(1,:) = cmap(1,:); colorcode(numberOfConditions,:)=cmap(64,:);
    cmapstep = 64/(numberOfConditions-1);
    for i = 1:numberOfConditions-2
        colorcode(i+1,:)=cmap(round(cmapstep*i),:);
    end
end