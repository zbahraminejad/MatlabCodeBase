%% Name_change_IX_Micro_1site_3color_timelapse.m
% This script converts the filenames of image files generated using ix
% microXL into a more useful format
% The script assumes that image files are stored within subfolders named
% 'Timepoint_XX' of the folder 'dir'. The files are renamed and moved into
% 'dir'. The output format is 'Row_Column_Site_Channel_Timepoint.tif'
% Written by Steve Cappell
% 120919

%% Init
clear
close all
clc

time1=tic;  %Starts Timer
%% Load Directory with images
dir='/Volumes/CAPPELL-2TB/20130323-MCF10A-EGF/1121';     %'/Users/scappell/Documents/Meyer_Lab/Data/2012-07-04/32';
alldir=getSubdirectories(dir); %Creates list of all the subfolders

%% Designate Channels
channels={'CFP','YFP','RFP'};  %Change this to match the correct channels used in the experiment  %e.g. {'Channel1','Channel2','Channel3'}  

for i=1:length(alldir);
    disp(i)                                                                               %Displays what timepoint you are on
    time2=tic;                                                                            %Starts timer to record time of each timepoint
    filenames=getFilenames([dir,'/TimePoint_',num2str(i)]);                           %Returns list of all file names inside the folder
    filenames=filenames(boolRegExp(filenames,'\.tif$') & ~boolRegExp(filenames,'thumb')); %removes accessory files and thumbnails from the file list
    
    %% extract wellname, site, and channel from original file name
    for j=1:length(filenames);
        wellname=char(getTokens(filenames(j),'^[^_]+_([A-Z][0-9][0-9])_'));  %Finds well name in format wellColumn (eg B04)
        row_column=wellNameToRowColumn(wellname);                            %Converts well name to format row_column (eg 2_4)
        site=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_s([0-9])_')); %Finds site image was taken (eg 1)
        channel=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_s[0-9]_w([1-9])'));  %Finds channel number (eg 1)
        
        
        %% move file
         newFileName=[row_column,'_',site,'_',channels{str2double(channel)},'_',num2str(i),'.tif'];  %designates new name (eg 2_4_1_EYFP_1)
%          newFileName=[wellname,'_1_',channels{str2double(channel)},'_',num2str(i),'.tif'];  %designates new name (eg B04_1_EYFP_1) %Use this line instead of previous if you want a different file name format
         oldFileName=['TimePoint_',num2str(i),'/',filenames{j}];        %Designates old file name to be used in old file path
         movefile([dir,'/',oldFileName],[dir,'/',newFileName],'f'); %Moves file from old folder to new folder and changes the name
            
    end
    toc(time2) %Displays time it takes to rename all the files in each Time point
end
toc(time1)  %Displays total time to rename all files
