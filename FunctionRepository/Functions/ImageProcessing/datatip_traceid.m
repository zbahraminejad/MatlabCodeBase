function output_txt = datatip_traceid(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

global traceid

%%%%% get x and y coordinates for frame & signal intensity
pos = get(event_obj,'Position');
frame = pos(1);
signal = num2str(pos(2),4);
output_txt = {['X: ',num2str(frame)],['Y: ',signal]};

%%%%% get x and y coordinates for frame & signal intensity
tar = get(event_obj,'Target');
traceid = get(tar,'DisplayName');  %returns cellid #