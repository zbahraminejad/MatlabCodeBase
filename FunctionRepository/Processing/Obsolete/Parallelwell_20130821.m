format short;
cd('H:\Documents\Timelapse\Code\Development\Processing')       % change to directory where program resides

time1=tic;

%rows = {'B' 'C' 'D' 'E' 'F' 'G' 'H'}; %alphabetical character
%cols = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};  %two numeric characters

%%% 20130715
%rows = {'B' 'C' 'D' 'E'};
%cols = {'03' '04'};
%sites = {'1' '2' '3' '4'};

%%% 20130719
%rows = {'B' 'C' 'D' 'E'};
%cols = {'04' '05' '06' '07' '08' '09' '10'};
%sites = {'1' '2' '3' '4'};

rows = {'B' 'C' 'D' 'E'};
cols = {'06' '08' '09'};
sites = {'1' '2' '3' '4'};

manualwells = {
    'B' '03' '1';
    'B' '03' '2';
    'B' '03' '3';
    'B' '03' '4';
    'B' '04' '1';
    'B' '04' '2';
    'B' '04' '4';
    };

manualcontrol = 0;
numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
shots=numrows*numcols*numsites;
if manualcontrol==1
    shots=size(manualwells,1);
end
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
    
    %%% Timelapse pipeline 08/15/2013 %%%%%%%%%%%%%%%%%%%%%%%%
    %Timelapse_1_formatfiles(row,col,site)
    %Timelapse_2_extractfeatures_test_postbs(row,col,site);
    %<no parallel> Timelapse_3_detectbadframes;
    %Timelapse_4_calculatejitters(row,col,site);
    %%%%% Immunostain pipeline can begin now %%%%%
    %Timelapse_5_trackcells(row,col,site);
    %Timelapse_6_gatetraces(row,col,site);
    
    %%% Immunostain pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %<no parallel> Immunostain_1_formatfiles;
    %<no parallel> Immunostain_2_calcbleedthroughrate(row,col,site); %only if extrapolating
    %Immunostain_3_updatewellsssjitter(row,col,site);
    %%% continue on to Timelapse_5_trackcells %%%%%%%
    
    %%% Timelapse pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Timelapse_0_formatfiles(row,col,site)
    %Timelapse_1_extractfeatures(row,col,site);
    %Timelapse_2_trackcells(row,col,site);
    %Timelapse_3_gatetraces(row,col,site);
end
toc(time1)
cd ..

%%% previous pipeline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A_format(row,col,site);
%A_crop_directory(row,col,site);
%B_extractfeatures_directory(row,col,site);
%B_extractfeatures_directory_nocdt1(row,col,site);
%C_trackcells_directory(row,col,site);
%D_tracesignals(row,col,site);