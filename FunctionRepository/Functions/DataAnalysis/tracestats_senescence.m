function [goodsenescence,badsenescence] = tracestats_senescence(row,col,site,path)
%%% Purpose: Calculate percentage of cells that are never divide %%%

%row=[3];col=[3 4];site=1:4;
%path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';

[signal1,signal2,mitoses,~] = combinedata_realbirths_history(row,col,site,path);

tracelist = 1:size(signal1,1);
totalcellnum = length(tracelist);
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

cd ..\Analysis; %return to this directory