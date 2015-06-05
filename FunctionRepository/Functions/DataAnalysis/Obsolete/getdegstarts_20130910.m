function degstarts=getdegstarts(sampletraces)
samplesize=size(sampletraces,1);
degstarts=zeros(samplesize,1);
for i=1:samplesize
    signal=sampletraces(i,:);
    sig_revslope=getslope_reverse(signal,1:5);
    sig_height=-signal;
    sig_time=(1:length(signal))/length(signal);
    filter=smooth(sig_revslope+sig_height-sig_time);
    lowerleft=find(filter==max(filter));
    sig_shortrevslope=getslope_reverse(signal,1);
    diffssrv=[diff(sig_shortrevslope),0];
    degstarts(i)=find(diffssrv>0 & [1:length(signal)]<lowerleft,1,'last');
end
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
x=1:length(signal);
plot(x,signal);
hold on;
plot(x,filter,'r');
%}