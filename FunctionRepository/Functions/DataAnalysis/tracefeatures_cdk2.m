%function [mitosistime,g1entry,g0toolong] = tracefeatures_cdk2(row,col,measurelimit,drugspike,path)
row=[3];col=[3 4];site=1:4;
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
measurelimit=20*5;drugspike=0;
%%%% Plot signal after post-drugspike mitosis vs drug-mitosis delay %%%%%%%%%%%%%%%%%
% USAGE: - row and col: specify well
%        - measurelimit: max length (in frames) for g1 or g1 latency
%        - drugspike: time (in frames) of drug addition
%        - featurechoice: 1 = G1 length; 2 = G1 latency
%%%%%%%%%%%%%% viewer settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1min = 0; y1max = 2; y1step = 0.4;
y2min = 0; y2max = 2; y2step = 0.4;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses,firstunique] = combinedata_realbirths_history(row,col,site,path);

%%%% gate out traces of interest %%%%%%%%%%%%%%%
tracelist = 1:size(signal1,1);
totalcellnum = length(tracelist);
numframes = size(signal1,2);
badsensorscreen = find(sum(signal2>0.1,2)<20 | min(signal2,[],2)>0.01);     %specifically for PIP degron sensor
numbad = length(badsensorscreen);
if numbad==totalcellnum    %assume this is the case only for cell lines w/ no sensor
    badsensorscreen = [];
    numbad = 0;
end
numgood = length(tracelist)-numbad;
nomitosesscreen = find(max(mitoses,[],2)==0);
numsenescent = length(nomitosesscreen);
badsensornomito = sum(ismember(badsensorscreen,nomitosesscreen));
goodsensornomito = numsenescent-badsensornomito;
badsenescence = round(100*badsensornomito/numbad);
goodsenescence = round(100*goodsensornomito/numgood);
%fprintf('bad sensor senescence = %0.0f  (n = %0.0f)\n',badsenescence,numbad);
%fprintf('good sensor senescence = %0.0f  (n = %0.0f)\n',goodsenescence,numgood);
senescence = goodsenescence;
screenlist = unique([badsensorscreen;nomitosesscreen]);
tracelist(ismember(tracelist,screenlist)) = [];
cdk2 = signal1(tracelist,:);
sensor2 = signal2(tracelist,:);
samplesize = length(tracelist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

firstunique = firstunique(tracelist,:);
xtime = 1:numframes;
measuremargin = 10;   %this is the margin needed to identify cdt1 peak or hdhb rise
maxmitosis = numframes-measurelimit-measuremargin;  %latest possible mitosis
goodmitoses = mitoses(tracelist,:);
tempzeros = zeros(samplesize,10);
mitosesafterspike = tempzeros;    %1st mitosis after drugspike

IMT = tempzeros;
IMTtoolong = tempzeros;

sentry = tempzeros;
sentryactual = tempzeros;
Sflag = tempzeros;
g1toolong = tempzeros;

g1entry = tempzeros;
g1entryactual = tempzeros;
G1flag = tempzeros;
g0toolong = tempzeros;

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
        %%% find next mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if mitosesafterspike(i,k+1)
            IMT(i,k) = mitosesafterspike(i,k+1)-curpoint;
            nextpoint = mitosesafterspike(i,k+1)-30; %for finding g1entry point, shortest G1/S to M is 8hrs
        else
            nextpoint = numframes;
            IMTtoolong(i,k) = 1;
        end
        %%% find approximate S-phase entry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        smoothcdk2 = smooth(cdk2(i,:));
        window = smoothcdk2(curpoint+1:nextpoint);
        sentry(i,k) = find(window>1.0,1,'first');
        if isempty(sentry(i,k))
            g1toolong(i,k) = IMTtoolong(i,k)==1;  %G1 actually too long
            Sflag(i,k) = IMTtoolong(i,k)==0; %error in CDK2 trace
            sentry(i,k)=nextpoint-curpoint; %search for G1 entry through whole cycle
        end
        sentrytime = curpoint+sentry(i,k);
        sentryactual(i,k) = sentrytime;   %for visualization only
        %%% find G1 entry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        g1slope = getslope_forward(smoothcdk2',1:10);
        g1depth = 0.5-smoothcdk2';
        g1dist = 1:length(smoothcdk2);
        risefilter = smooth(g1slope*0.5+g1depth+g1dist*0.003);
        risefilter = risefilter(curpoint+1:curpoint+sentry(i,k));
        g1entry(i,k) = find(risefilter==max(risefilter),1,'first');
        g1entrytime = curpoint+g1entry(i,k);
        g1entryactual(i,k) = g1entrytime;   %for visualization only
        boundarycheck = g1entryactual(i,k)>=numframes-measuremargin;  %could only happen for g0toolong
        G1flag = IMTtoolong(i,k)==1 && smoothcdk2(end)<0.75;
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
noisiness=mean(abs(diff(cdk2,1,2)),2)*100;
leavein(noisiness>noisethresh,:)=0;

%%%% subset toggling (uncomment to activate ) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leavein(g1toolong==1)=0;

%%%% label quiescent traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1toolong(sentry>measurelimit)=1;        %ignore risetime measurements
g0toolong(g1entry>measurelimit)=1;

%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
showsignal2 = 0;
if panelvisualize
    sample=find(sum(leavein,2)); %finds all traces (rows) that have valid data
    %sample=sample(find(sample>=192 & sample<=217));
    %sample=sample(find(sample>=180 & sample<=200));
    %for i=1:samplesize
    %sample=[8 20 34 35 52 54 56 65 98 110 183 184 188 192 197 200];
    for cc=1:length(sample)
        i=sample(cc);
        figure(ceil(cc/20));             %only 24 plots per figure
        set(gcf,'color','w');
        set(gcf,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]);
        %subaxis(4,5,mod(cc-1,20)+1,'ML',0.02,'MR',0.02,'MT',0.03,'MB',0.03,'SH',0.03); %5x4 screen
        subaxis(4,5,mod(cc-1,20)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
        
        signal1trace = cdk2(i,:);
        signal2trace = sensor2(i,:);
        %%%% uncomment to visualize filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        smoothsignal = smooth(signal1trace);
        g1slope = getslope_forward(smoothsignal',[1:10]);
        g1depth = 0.5-smoothsignal';
        g1dist = 1:length(smoothsignal);
        signal2trace = smooth(g1slope*0.5+g1depth+g1dist*0.003);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if showsignal2==1
            [haxes,hline1,hline2] = plotyy(xtime,signal1trace,xtime,signal2trace);
            axes(haxes(1));
        else
            hline1=plot(xtime,signal1trace);
        end
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
        scatter(firstunique(i),signal1trace(firstunique(i)),120,'ko','markerfacecolor','k');
        scatter(sentryactual(i,markable),signal1trace(sentryactual(i,markable)),'go','markerfacecolor','g');
        scatter(g1entryactual(i,markable),signal1trace(g1entryactual(i,markable)),120,'ro','markerfacecolor','r');
        if showsignal2==1
            axes(haxes(2)); hold on;
            axis([0 numframes y2min y2max]);
            set(gca,'YTick',y2min:y2step:y2max);
            set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
            %line([xtime(1) xtime(end)],[0 0],'Color','k');
        end
    end
end

%%%%%%%%%%%%%%%%%% Compile all good data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mitosistime = mitosesafterspike(logical(leavein));
IMT = IMT(logical(leavein));
sentry = sentry(logical(leavein));
g1entry = g1entry(logical(leavein));
g1toolong = g1toolong(logical(leavein));
g0toolong = g0toolong(logical(leavein));
IMTtoolong = IMTtoolong(logical(leavein));

%set(gcf,'color','w','PaperPosition',[0 0 18 12]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);