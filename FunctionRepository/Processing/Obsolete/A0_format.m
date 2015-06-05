%% Name_change_IX_Micro_1site_3color_timelapse.m
% This script converts the filenames of image files generated using ix
% microXL into a more useful format
% The script assumes that image files are stored within subfolders named
% 'Timepoint_XX' of the folder 'dir'. The files are renamed and moved into
% 'dir'. The output format is 'Row_Column_Site_Channel_Timepoint.tif'
% Written by Steve Cappell
% 120919
cd ..\Functions; %change directory for function calls

time1=tic;  %Starts Timer
%% Load Directory with images
dir='H:\Documents\Timelapse\Timescape\20130514_MCF10A-p21KO_24hr_6min\';
preformatdir=[dir,'ix_micro'];
postformatdir=[dir,'Raw'];
alldir=getSubdirectories(preformatdir); %Creates list of all the subfolders

%% Designate Channels
channels={'CFP','YFP'};  %Change this to match the correct channels used in the experiment  %e.g. {'Channel1','Channel2','Channel3'}  

for i=1:length(alldir);
    disp(i)                                                                               %Displays what timepoint you are on
    time2=tic;                                                                            %Starts timer to record time of each timepoint
    filenames=getFilenames([preformatdir,'\TimePoint_',num2str(i)]);                           %Returns list of all file names inside the folder
    filenames=filenames(boolRegExp(filenames,'\.tif$') & ~boolRegExp(filenames,'thumb')); %removes accessory files and thumbnails from the file list
    
    %% extract wellname, site, and channel from original file name
    for j=1:length(filenames);
        wellname=char(getTokens(filenames(j),'^[^_]+_([A-Z][0-9][0-9])_'));  %Finds well name in format wellColumn (eg B04)
        wellname=[wellname(1),'_',wellname(2:3),'_1'];
        %row_column=wellNameToRowColumn(wellname);                            %Converts well name to format row_column (eg 2_4)
        %site=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_s([0-9])_')); %Finds site image was taken (eg 1)
        channel=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_w([1-9])'));  %Finds channel number (eg 1)      
       
        %% move file
         oldFileName=['TimePoint_',num2str(i),'\',filenames{j}];        %Designates old file names to be used in old file path
         newFileName=[wellname,'\',channels{str2double(channel)},'_',num2str(i),'.tif'];  %designates new name (eg 2_4_1_EYFP_1)
         if i==1
             system(['mkdir ',postformatdir,'\',wellname]);
         end
         movefile([preformatdir,'\',oldFileName],[postformatdir,'\',newFileName],'f');     %Moves files from old folder to new folder and changes the names
            
    end
    toc(time2) %Displays time it takes to rename all the files in each Time point
end
toc(time1)  %Displays total time to rename all files
