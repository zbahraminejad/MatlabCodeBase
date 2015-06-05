function [rQ, rM, labelmask] = regionQM(mask, image, lo, hi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mask, black/white image or labelled%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%output will have the same label as input, cell identity maintained%%%%%%%%
%logical mask is much faster than label mask%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image, gray scale image for analysis, same dimension as mask%%%%%%%%%%%%%%
%lo, low end of quantile, must be between [0, 1]%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hi, high end of quantile, must be between [0, 1]%%%%%%%%%%%%%%%%%%%%%%%%%%
%if lo==hi, rQ returns the lo quantile of each ROI%%%%%%%%%%%%%%%%%%%%%%%%%
%if lo<hi, rQ returns the mean intensity of each ROI between lo&hi quantile
%rM returns the mode of intensity for each ROI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%labelmask returns the labeled mask, same as bwlabel(mask)%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if max(mask(:)) == 1
        CC = bwconncomp(mask);
        imageSize = CC.ImageSize;
        numObjs = CC.NumObjects;
        rQ = zeros(numObjs, 1);
        rM = zeros(numObjs, 1);
        labelmask = zeros(imageSize);
    
        for k = 1:numObjs
            IdxList = cell2mat(CC.PixelIdxList(k));
            Intensity = image(IdxList);
            rM(k) = mode(Intensity);
            if lo == hi
                rQ(k) = quantile(Intensity, lo);
            else
                Ilo = quantile(Intensity, lo);
                Ihi = quantile(Intensity, hi);
                rQ(k) = mean(Intensity(Intensity>=Ilo & Intensity<=Ihi));
            end
            labelmask(IdxList) = k;
        end
        
    else
        imageSize = size(mask);
        Objs = sort(unique(mask(:)));
        numObjs = length(Objs)-1;
        rQ = zeros(numObjs, 1);
        rM = zeros(numObjs, 1);
        labelmask = zeros(imageSize);
    
        for k = 1:numObjs
            IdxList = find(mask == Objs(k+1));
            Intensity = image(IdxList);
            rM(k) = mode(Intensity);
            if lo == hi
                rQ(k) = quantile(Intensity, lo);
            else
                Ilo = quantile(Intensity, lo);
                Ihi = quantile(Intensity, hi);
                rQ(k) = mean(Intensity(Intensity>=Ilo & Intensity<=Ihi));
            end
            labelmask(IdxList) = Objs(k+1);
            % A = abs(labelmask - mask);
            % sum(A(:)) shall equal to zero. if not, idxlist is wrong
        end
 
    end
end
