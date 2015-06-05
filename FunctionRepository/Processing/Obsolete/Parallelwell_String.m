%%% template
%rows = {'B' 'C' 'D' 'E' 'F' 'G' 'H'}; %alphabetical character
%cols = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};  %two numeric characters

%%% 20130715
% rows = {'B' 'C' 'D' 'E'};
% cols = {'03' '04' '05' '06'};
% sites = {'1' '2' '3' '4'};

%%% 20130719
%rows = {'B' 'C' 'D' 'E'};
%cols = {'04' '05' '06' '07' '08' '09' '10'};
%sites = {'1' '2' '3' '4'};

%%% 20130831
% rows={'B' 'C' 'D' 'E' 'F' 'G'};
% cols={'02' '03' '04' '05' '06' '07' '08' '09' '10'};
% sites={'1' '2'};

%%% 20130905
% rows={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
% cols={'02' '03' '04' '05' '06' '07' '08' '09' '10'};
% sites={'1' '2'};

%%% 20130920
% rows = {'B' 'C' 'D' 'E' 'F' 'G'};
% cols = {'02' '03' '04' '05' '06' '07' '08' '09'};
% sites = {'1' '2'};

%%% JY CycD1-FKBP
% rows = {'B'};
% cols = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
% sites = {'1' '2'};

%%% SS 20121211 MCF10A DHFRp21 & MLN
% rows = {'4' '5'};
% cols = {'3' '5'};
% sites = {'0' '1'};

%%% SS 20130222 MCF10Ap21null DHFRp21
% rows = {'6' '8'};
% cols = {'5' '6' '7' '8'};
% sites = {'1' '2'};

%%% 20131116 CDK4_pRb_Abs S780 MvsG
rows = {'1' '2' '3' '4' '5' '6' '7' '8'};
cols = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
sites = {'1' '2' '3' '4'};

%%% sample
% rows = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
% cols = {'08'};
% sites = {'1' '2'};

manualwells = {
    'B' '04' '1'
    'B' '06' '1'
    'B' '08' '1'
    'C' '04' '1'
    'C' '06' '1'
    'C' '08' '1'
    };

manualcontrol = 0;
numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
shots=numrows*numcols*numsites;
if manualcontrol==1
    shots=size(manualwells,1);
end
time1=tic;
parfor shot=1:shots
    if manualcontrol==1
        row=manualwells{shot,1};
        col=manualwells{shot,2};
        site=manualwells{shot,3};
    else
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites{siteidx};
        colidx=mod(ceil(shot/numsites),numcols);
        if colidx==0
            colidx=numcols;
        end
        col=cols{colidx};
        rowidx=ceil(shot/(numcols*numsites));
        row=rows{rowidx};
    end
    fprintf([row,'_',col,'_',site,'\n']);
    
    %%% Format Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FormatFiles(row,col,site);
    %FormatFiles_Rename(row,col,site);
    %FormatFiles_Immunostain(row,col,site);
    
    %%% Timelapse %%%%%%%%%%%%%%%%%%%%%%%%%
    %Timelapse_merge(row,col,site);
    
    %%% Immunostain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Immunostain_1_CalcBleedthroughRate(row,col,site);
    %Immunostain_2_AddToTimelapse(row,col,site);
    
end
toc(time1)