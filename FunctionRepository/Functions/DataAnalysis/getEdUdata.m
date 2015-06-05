%function [sensorcalc,EdU] = getEdUdata(row,col,site)
row=[4 5];col=7;site=1:4;
%%%% Plot signal after post-drugspike mitosis vs drug-mitosis delay %%%%%%%%%%%%%%%%%
% USAGE: - row and col: specify well
%        - measurelimit: max length (in frames) for g1 or g1 latency
%        - drugspike: time (in frames) of drug addition
%        - featurechoice: 1 = G1 length; 2 = G1 latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Directory Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses,firstunique,EdU,DAPI] = combinedata_EdU(row,col,site,path);

%%%% gate out traces of interest %%%%%%%%%%%%%%%
tracelist = 1:size(signal1,1);
numframes = size(signal1,2);
screen = find(sum(signal2>0.1,2)<20 | min(signal2,[],2)>0.01);
tracelist(ismember(tracelist,screen)) = [];
sensor = signal2(tracelist,:);
EdU = EdU(tracelist);
dapi = DAPI(tracelist);
samplesize = length(tracelist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% viewer settings %%%%%%%%%%%%%%%%%
y1min = 0; y1max = 2.5; y1step = 0.5;
y2min = 0; y2max = 0.8; y2step = 0.2;
signalwidth = 2;
framesperhr = 5;
xtime=1:numframes;
realtime=xtime/framesperhr;

%%%% calculate final sensor value %%%%%%%%%%%%%%
%sensorcalc = mean(sensor(:,end-2:end-1),2);  %avg of last 15 minutes (excluding stain)
sensorcalc = sensor(:,end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sensor value vs EdU incorp %%%%%%%%%%%%%%%%%
%{
scatter(EdU,sensorcalc,40);xlabel('Cyclin D1 Immunofluorescence');ylabel('Sensor Level');
scatter(sensorcalc,EdU,40);ylabel('Cyclin D1 Immunofluorescence');xlabel('Sensor Level');
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}

%%%% for low vals, find most recent fall %%%%%%%%%%%%%%%%%%%
lowset = find(sensorcalc<0.05);
lowpoi = zeros(samplesize,1);
lowfilter = zeros(samplesize,numframes);
for i=lowset'
    signal = sensor(i,:);
    reversesignal = signal(length(signal):-1:1);
    rs_slope = getslope_forward(reversesignal,1:10);
    rs_height = 0.5-reversesignal;
    rs_time = 1:length(signal);
    reversefilter = smooth(rs_slope*0.25+rs_height-rs_time*0.001);
    lowfilter(i,:) = reversefilter(length(reversefilter):-1:1);
    lowpoi(i) = find(lowfilter(i,:)==max(lowfilter(i,:)),1,'last');
end
timesincedrop = numframes-lowpoi;
%{
scatter(timesincedrop(lowset)/5,EdU(lowset),40);
xlabel('Time since signal drop (hr)');
ylabel('Cyclin D1 Immunofluorescence');
axis([0 15 0 3.5]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
%{
xmin=min(timesincedrop(lowset)); xmax=max(timesincedrop(lowset));
x=xmin:xmax;
y=zeros(1,length(x));
e=zeros(1,length(x));
for i=1:length(x)
    index=find(timesincedrop(lowset)==x(i));
    y(i)=mean(EdU(lowset(index)));
    e(i)=std(EdU(lowset(index)));
end
x=x/5;
errorbar(x,y,e,'o','markeredgecolor','k','markerfacecolor',[.49 1 .63],'linewidth',1.5);
xlabel('Time since signal drop (hr)');
ylabel('EdU Incorporation');
axis([0 15 0 3.5]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}

%%%% for high vals, find most recent inflection %%%%%%%%%%%%%%%%%%%
highset = find(sensorcalc>0.1);
highpoi = zeros(samplesize,1);
for i=highset'
    highpoi(i) = find(sensor(i,:)<0.05,1,'last');
end
timesincerise = numframes-highpoi;
%{
scatter(timesincerise(highset)/5,EdU(highset),40,'markeredgecolor','k','markerfacecolor',[.49 1 .63],'linewidth',1.5);
scatter(timesincerise(highset)/5,EdU(highset),40);
xlabel('Time since signal rise (hr)');
ylabel('EdU Incorporation');
axis([0 45 0 3.5]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}

%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 0;
showsignal2 = 0;
if panelvisualize
    %sample=1:samplesize;
    sample=highset;
    for cc=1:length(sample)
        i=sample(cc);
        figure(ceil(cc/20));             %only 24 plots per figure
        set(gcf,'color','w');
        set(gcf,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]);
        subaxis(4,5,mod(cc-1,20)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
        sensortrace = sensor(i,:);
        if showsignal2==1
            signal2trace = lowfilter(i,:);
            [haxes,hline1,hline2] = plotyy(xtime,sensortrace,xtime,signal2trace);
            axes(haxes(1));
        else
            hline1=plot(xtime,sensortrace);
        end
        hold on;
        axis([0 numframes y1min y1max]);
        set(gca,'YTick',y1min:y1step:y1max);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        
        %scatter(numframes,EdU(i),120,'ro','markerfacecolor','r');
        %scatter(lowpoi(i),sensortrace(lowpoi(i)),80,'ro','markerfacecolor','r');
        scatter(highpoi(i),sensortrace(highpoi(i)),80,'ro','markerfacecolor','r');
        
        if showsignal2==1
            axes(haxes(2)); hold on;
            axis([0 numframes y2min y2max]);
            set(gca,'YTick',y2min:y2step:y2max);
            set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
            %line([xtime(1) xtime(end)],[0 0],'Color','k');
            %}
        end
    end
end


%{
scatter(EdU,sensorcalc);
xlabel('EdU Incorporation');
ylabel('p21-sensor');
axis([-0.5 3.5 -0.1 2]);
set(gcf,'color','w','PaperPosition',[0 0 6 9]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}

cd ..\Analysis; %return to this directory