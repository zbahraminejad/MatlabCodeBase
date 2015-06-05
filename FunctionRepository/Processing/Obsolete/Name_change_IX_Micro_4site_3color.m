%% Init
clear
close all
clc

time1=tic;  %Starts Timer
%% Load Directory with images
dir='L:\External JY\JY\JY-BJ5-CycDFKBP-DHFRp21-061213\2013-06-12\1512\';     %'/Users/scappell/Documents/Meyer_Lab/Data/2012-07-04/32';
newdir='L:\JY-BJ5-CycDFKBP-DHFRp21-061213';
mkdir(newdir)
alldir=getSubdirectories(dir); %Creates list of all the subfolders

%% Designate Channels
channels={'ECFP','Texas Red','YFP'};  %Change this to match the correct channels used in the experiment  %e.g. {'Channel1','Channel2','Channel3'}  

for i=1:length(alldir);
    disp(i)                                                                               %Displays what timepoint you are on
    time2=tic;                                                                            %Starts timer to record time of each timepoint
    filenames=getFilenames([dir,'\TimePoint_',num2str(i)]);                           %Returns list of all file names inside the folder
    filenames=filenames(boolRegExp(filenames,'\.tif$') & ~boolRegExp(filenames,'thumb')); %removes accessory files and thumbnails from the file list
    
    %% extract wellname, site, and channel from original file name
    for j=1:length(filenames);
        wellname=char(getTokens(filenames(j),'^[^_]+_([A-Z][0-9][0-9])_'));  %Finds well name in format wellColumn (eg B04)        
        site=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_s(.)_')); %Finds site image was taken (eg 1)
        channel=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_s._w([1-9])'));  %Finds channel number (eg 1)      
        %% create new folder for well and site
        if ~exist([newdir,'\', wellname],'dir')
            mkdir([newdir,'\', wellname]);
        end
        
        if ~exist([newdir,'\', wellname,'\site' num2str(site)],'dir')
            mkdir([newdir,'\', wellname,'\site' num2str(site)]);
        end
        
        %% move file
         newFileName=[wellname,'_site',site,'_',channels{str2double(channel)},'_',num2str(i),'.tif'];  %designates new name (eg 2_4_1_EYFP_1)
         oldFileName=['Timepoint_',num2str(i),'\',filenames{j}];        %Designates old file name to be used in old file path
         movefile([dir,'\',oldFileName],[newdir,'\', wellname,'\site' num2str(site),'\',newFileName],'f'); %Moves file from old folder to new folder and changes the name
            
    end
    toc(time2) %Displays time it takes to rename all the files in each Time point
end
toc(time1)  %Displays total time to rename all files
