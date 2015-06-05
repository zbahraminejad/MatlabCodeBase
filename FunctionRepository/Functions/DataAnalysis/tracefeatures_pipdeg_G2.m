%function [g2entry,mitosistime,g2toolong] = getcyclelengths_pipdegron_G2(row,col,measurelimit,drugspike)
row=[2 3 4];col=[3 4];site=1:4;
measurelimit=20*5;drugspike=0;
%%%% Plot signal after post-drugspike mitosis vs drug-mitosis delay %%%%%%%%%%%%%%%%%
% USAGE: - row and col: specify well
%        - measurelimit: max length (in frames) for g1 or g1 latency
%        - drugspike: time (in frames) of drug addition
%        - featurechoice: 1 = G1 length; 2 = G1 latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Directory Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
%%%%%%%%%%%%%% viewer settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1min = 0; y1max = 2.5; y1step = 0.5;
%y2min = 0; y2max = 1; y2step = 0.2;
y2min = 0; y2max = 0.8; y2step = 0.2;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses,firstunique] = combinedata_realbirths_history(row,col,site,path);

%%%% gate out traces of interest %%%%%%%%%%%%%%%
tracelist = 1:size(signal1,1);
numframes = size(signal1,2);
badsensorscreen = find(sum(signal2>0.1,2)<20 | min(signal2,[],2)>0.01);
nomitosesscreen = find(max(mitoses,[],2)==0);
screenlist = unique([badsensorscreen;nomitosesscreen]);
tracelist(ismember(tracelist,screenlist)) = [];
angie = signal1(tracelist,:);
sensor = signal2(tracelist,:);
samplesize = length(tracelist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstunique = firstunique(tracelist,:);
xtime = 1:numframes;
measuremargin = 10;   %this is the margin needed to identify cdt1 peak or hdhb rise
maxmitosis = numframes-measurelimit-measuremargin;  %latest possible mitosis
goodmitoses = mitoses(tracelist,:);
tempzeros = zeros(samplesize,10);
mitosesafterspike = tempzeros;    %1st mitosis after drugspike
Sflag = tempzeros;
G1flag = tempzeros;
g1toolong = tempzeros;
g0toolong = tempzeros;
pipdegdrop = tempzeros;
pipdegdropactual = tempzeros;
threshpercent = 0.15;
g1entry = tempzeros;
g1entryactual = tempzeros;
fallfilter = zeros(size(angie));
g0entry = tempzeros;
g0entryactual = tempzeros;
leavein = ~tempzeros;   %ones
maxprev = 20;
set(0,'Units','pixels');           %sets screensize units by pixels
screendims = get(0,'ScreenSize');  %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);

for i=1:samplesize
    temp = goodmitoses(i,:);
    temp = sort(temp(temp>10));
    mitosesafterspike(i,1:length(temp)) = temp;
end

for i=1:samplesize
    if mitosesafterspike(i,1)==0
        continue        %skip to next iteration.  any traces without postspike mitosis ignored
    end
    k = 1;
    while mitosesafterspike(i,k)
        curpoint = mitosesafterspike(i,k);
        if mitosesafterspike(i,k+1)
            nextpoint = mitosesafterspike(i,k+1);
        else
            nextpoint = numframes;
        end
        smoothsensor = smooth(sensor(i,:));
        pipdegg1 = smoothsensor(curpoint+1:nextpoint);
        [~,pipdegdrop(i,k)] = find(pipdegg1<0.05,1,'first');
        pipdegdroptime = curpoint+pipdegdrop(i,k);
        pipdegdropactual(i,k) = pipdegdroptime;    %only for visualization
        lowercheck = pipdegdroptime<=curpoint+2;
        uppercheck = pipdegdroptime>=nextpoint-measuremargin;
        boundarycheck = lowercheck || uppercheck;
        g1toolong(i,k) = boundarycheck && nextpoint==numframes;
        Sflag(i,k) = boundarycheck && nextpoint<numframes;
        
        %%% find G1 Entry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        smoothangie = smooth(angie(i,:));
        g1slope = getslope_forward(smoothangie',[1:10]);
        g1depth = 0.5-smoothangie';
        g1dist = 1:length(smoothangie);
        risefilter = smooth(g1slope*0.5+g1depth+g1dist*0.003);
        if g1toolong(i,k)==1
            risefilter = smooth(g1slope+g1depth+g1dist*0.003);
            pipdegdroptime = numframes;
        end
        risefilter = risefilter(curpoint+1:pipdegdroptime);
        g1entry(i,k) = find(risefilter==max(risefilter),1,'first');
        g1entrytime = curpoint+g1entry(i,k);
        g1entryactual(i,k) = g1entrytime;   %for visualization only
        boundarycheck = g1entryactual(i,k)>=numframes-measuremargin;  %could only happen for g1toolong
        G1flag = g1toolong(i,k)==1 && smoothangie(end)<0.75;
        g0toolong(i,k) = boundarycheck || G1flag;
        
        k = k+1;
    end
end

%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%%%%%%%%%%%%%%%%%%
leavein(mitosesafterspike==0)=0;
leavein(mitosesafterspike>=maxmitosis)=0;
leavein(logical(Sflag))=0;  %errant trace
leavein(g1entryactual<drugspike)=0;

%%%% gate out noisy hDHB traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisethresh = 30;
noisiness=mean(abs(diff(angie,1,2)),2)*100;
leavein(noisiness>noisethresh,:)=0;

%%%% subset toggling (uncomment to activate ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leavein(g1toolong==1)=0;

%%%% label quiescent traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1toolong(pipdegdrop>measurelimit)=1;        %ignore risetime measurements

%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
if panelvisualize
    sample=find(sum(leavein,2));
    %sample=sample(find(sample>=192 & sample<=217));
    sample=sample(find(sample>=180 & sample<=200));
    %for i=1:samplesize
    %sample=[8 20 34 35 52 54 56 65 98 110 183 184 188 192 197 200];
    for cc=1:length(sample)
        i=sample(cc);
        figure(ceil(cc/20));             %only 24 plots per figure
        set(gcf,'color','w');
        set(gcf,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]);
        %subaxis(4,5,mod(cc-1,20)+1,'ML',0.02,'MR',0.02,'MT',0.03,'MB',0.03,'SH',0.03); %5x4 screen
        subaxis(4,5,mod(cc-1,20)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
        signal1trace = angie(i,:);
        %signal2trace = 0;
        %signal2trace = smooth(cdt1(i,:));
        %signal2trace = fallfilter(i,:);
        
        %%%% uncomment to visualize filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        smoothsignal = smooth(signal1trace);
        g1curve = getcurvature_array(smoothsignal',[1:10]);
        g1depth = 1-smoothsignal';
        g1dist = 1:length(smoothsignal);
        signal2trace = smooth(g1curve+g1depth-g1dist*0.003);
        %}
        
        smoothsignal = smooth(signal1trace);
        g1slope = getslope_forward(smoothsignal',[1:10]);
        g1depth = 0.5-smoothsignal';
        g1dist = 1:length(smoothsignal);
        signal2trace = smooth(g1slope*0.5+g1depth+g1dist*0.003);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [haxes,hline1,hline2] = plotyy(xtime,signal1trace,xtime,signal2trace);
        axes(haxes(1)); 
        %hline1=plot(xtime,signal1trace);
        hold on;
        axis([0 numframes y1min y1max]);
        set(gca,'YTick',y1min:y1step:y1max);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        
        if drugspike
            line([drugspike drugspike],[y1min y1max],'Color','r');  %add vertical line when drug is added
        end
        
        markable = find(leavein(i,:));
        scatter(mitosesafterspike(i,markable),signal1trace(mitosesafterspike(i,markable)),120,'co','markerfacecolor','c');
        %scatter(g0entryactual(i,markable),signal1trace(g0entryactual(i,markable)),'go','markerfacecolor','g');
        scatter(g1entryactual(i,markable),signal1trace(g1entryactual(i,markable)),120,'ro','markerfacecolor','r');

        scatter(firstunique(i),signal1trace(firstunique(i)),120,'ko','markerfacecolor','k');
        
        
        
        axes(haxes(2)); hold on;
        axis([0 numframes y2min y2max]);
        set(gca,'YTick',y2min:y2step:y2max);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        %line([xtime(1) xtime(end)],[0 0],'Color','k');
        %}
        scatter(pipdegdropactual(i,markable),signal2trace(pipdegdropactual(i,markable)),120,'go','markerfacecolor','k');
    end
end

%%%%%%%%%%%%%%%%%% Compile all good data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%g1riselength = sentry-g1entry;        %length of g1 rise
mitosistime = mitosesafterspike(logical(leavein));
pipdegdrop = pipdegdrop(logical(leavein));       %goodg1length = g1length.*countg1;
g1entry = g1entry(logical(leavein));
%g1riselength = g1riselength(logical(leavein));
g1toolong = g1toolong(logical(leavein));
%quiescentcount = sum(sum(quiescent));

%set(gcf,'color','w','PaperPosition',[0 0 18 12]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);

cd ..\Analysis; %return to this directory