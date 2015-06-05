cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
nucr=16; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%%% setup tempwell %%%%%%%%%%%%%%%%%%
row=2;
col=1;
site=1;
track=779;
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
shottrack=[shot,'_track',num2str(track)];
load([datadir,shot,'_alldata'],'bestsp','best_rc','corruptlist','leaveoutlist');
load([datadir,shottrack,'attributes.mat'],'mon','msn','sbn','ssn','aor','anr','sor','snr','areabignuc','areasmallnuc');
framesperhr=5;
start=best_rc(track,1);
finish=best_rc(track,3);
finish=200;
frames=(start:finish);
time=frames/framesperhr; % calculate time (hr) instead of frame
time=frames;

%%%% set attribute
cytoconc = aor(start:finish);
nucconc= mon(start:finish);
attribute=cytoconc./nucconc;
attribute=nucconc;

%%%% Time Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'color','w')
plot(time,attribute,'color','r','linewidth',2);
%xlabel('Time (hrs)','fontsize',12);
xlabel('Frames','fontsize',12);
ylabel('signal','fontsize',12,'fontweight','b');
xlim([time(1) time(end)]);
%xlim([time(1) 40]);
ylim([min(attribute) max(attribute)]);
%ylim([500 1800]);
%ylim([0.6 1.6]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%  FFT calc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = fft(attribute(start:finish));
F(1)=[];
L=length(F);
nyquist = 1/2;
freq = (1:L/2)/(L/2)*nyquist*framesperhr;
P = abs(F(1:floor(L/2))).^2;
%%%% FFT Plot
subplot(2,1,2);
plot(freq,P,'linewidth',1);
title('Periodogram')
xlabel('Frequency (per hour)','fontsize',12);
ylabel('Power','fontsize',12,'fontweight','b');
xlim([0 2.5]);
ylim([0 2]);
noise=sum(P(8:38))     % anything above 2 hour cycles
noise=sum(P(4:38))     % anything above 4 hour cycles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

saveas(gcf,'h:\Downloads\Fig.jpg');
cd ..\Analysis; %return to this directory