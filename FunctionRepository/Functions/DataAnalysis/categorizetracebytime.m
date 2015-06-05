function [trace2high,trace2low]=categorizetracebytime(traces2,daughterstats,framereltomitosis,prctilethresh)
numtraces=size(traces2,1);
val2=ones(numtraces,1)*NaN;
for i=1:numtraces
    val2(i)=traces2(i,daughterstats(i,1)+framereltomitosis);
end
%figure,hist(val2,100);
trace2high=val2>prctile(val2,prctilethresh(1));
trace2low=val2<prctile(val2,prctilethresh(2));
end