function H = plotlinefill(varargin)

switch nargin
    case{1}
        signal = varargin{1};
        time = 1:size(signal,2);
        prop = 'r-';
        new_figure = true;
    case{2}
        time = varargin{1};
        signal = varargin{2};
        prop = 'r-';
        new_figure = true;
    case{3}
        time = varargin{1};
        signal = varargin{2};
        prop = varargin{3};
        new_figure = true;
    case{4}
        time = varargin{1};
        signal = varargin{2};
        prop = varargin{3};
        new_figure = varargin{4};
end

if ischar(new_figure)
    temp7 = regexpi(new_figure,'true');
    if isempty(temp7)
        new_figure = false;
        H = get(0,'CurrentFigure');
    else
        new_figure = true;
    end
elseif isnumeric(new_figure)
    H = new_figure;
    new_figure = false;
end

% prepare vectors for plotting
time2 = [time fliplr(time)];

sigm = nanmean(signal);
sigs =  nanstd(signal);
sigms = [sigm-sigs fliplr(sigm+sigs)];

sigmin = nanmin(signal);
sigmax = nanmax(signal);
sigmm = [sigmin fliplr(sigmax)];

%  check color
if strcmpi(prop(1),'r')
    c1 = [1 0 0];
elseif strcmpi(prop(1),'g')
    c1 = [0 1 0];
elseif strcmpi(prop(1),'b')
    c1 = [0 0 1];
else
    c1 = [1 1 1];
end

color1 = c1 + .7*(1-c1);
color2 = c1 + .8*(1-c1);

%
if length(prop) == 1
    prop(2) = '-';
end
if new_figure
    H = figure;
else
    figure(H);
end

whole_screen = get(0,'ScreenSize');
% max figure size - add
fig_size = whole_screen + [-4 -4+2*32 +4+4 4+4-2*32];
set(H,'OuterPosition',fig_size);

plot(time,sigm,prop,'LineWidth',1.5)
hold all
fill(time2(~isnan(sigms)), sigms(~isnan(sigms)),color1,'EdgeColor',color1,'FaceAlpha', 0.4);
fill(time2(~isnan(sigmm)), sigmm(~isnan(sigmm)),color2,'EdgeColor',color2,'FaceAlpha', 0.3);
legend([{'Mean'}  {'\pm Stddev'} {'Min/Max'} ],'Location','Best')