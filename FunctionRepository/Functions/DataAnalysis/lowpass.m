function filter_matrix = lowpass(order, cutoff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a lowpass filter that can be used with the filter2 function
% 
% Output: 
%     filter_matrix - matrix with filter coefficients that can be used with
%     the filter2 function
% Inputs:
%     order - The filter has size (order x order)
%	  cutoff - is the cutoff for the lowpass filter (cutoff specified between 0 and 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create desired frequency response
[f1,f2] = freqspace(order,'meshgrid');
d = find(f1.^2+f2.^2 < cutoff^2);
Hd = zeros(order);
Hd(d) = ones(size(d));

% Design the filter's impulse response
filter_matrix = fsamp2(Hd);