codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
datadir = [path,'Data\'];


%%% 20130715 %%%%%%%%%%%%%%%%%
%rows = {'B' 'C' 'D' 'E'};
%cols = {'03' '04' '05' '06'};
%sites = {'1' '2' '3' '4'};
%%% 20130719 %%%%%%%%%%%%%%%%%
%rows = {'B' 'C' 'D' 'E'};
%cols = {'04' '05' '06' '07' '08' '09' '10'};
%sites = {'1' '2' '3' '4'};

%rows = {'B' 'C' 'D' 'E'};
%cols = {'04' '05' '07' '10'};
%sites = {'1' '2' '3' '4'};
rows = {'D'}; cols = {'05'}; sites = {'4'};

fmax=230;
totalwells=length(rows)*length(cols)*length(sites);
emptymatrix=zeros(totalwells,fmax);
blurredmatrix=zeros(totalwells,fmax);
movieindex=cell(totalwells,1);
indexcount=0;
for rowidx = 1:length(rows)
    row = rows{rowidx};
    for colidx = 1:length(cols)
        col = cols{colidx};
        for siteidx = 1:length(sites)
            indexcount=indexcount+1;
            site = sites{siteidx};
            shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
            %fprintf('%s\n',shot);
            movieindex{indexcount}=shot;
            %%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            load([datadir,'wellsss_',shot,'.mat'],'wellsss');
            fmax=size(wellsss,3);
            %%% detect bad frames %%%%%%%%%%%%%%%%%%%%%%%%%
            numcells=zeros(1,fmax);
            mednucsize=zeros(1,fmax);
            for i=1:fmax
                numcells(i)=size(wellsss{i},1);
                mednucsize(i)=median(wellsss{i}(:,4));  %col 4 for nuc area
            end
            emptymatrix(indexcount,numcells==0)=1;
            %blurthresh=prctile(mednucsize,90)+6*iqr(mednucsize);
            blurthresh=median(mednucsize)*1.5;
            blurredmatrix(indexcount,mednucsize>blurthresh)=1;
        end
    end
end

%%% manually set missing frames for 20130715 experiment %%%
emptymatrix(strcmp(movieindex,'B_04_4'),[85 87 92 95 99 104 112 121 149 157])=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

badframes=emptymatrix | blurredmatrix;

fileID=fopen([datadir,'badframes.txt'],'w');
if sum(sum(emptymatrix))>0
    fprintf(fileID,'empty frames detected:\n');
    wellidx=find(sum(emptymatrix,2)>0);
    for i=wellidx'
        wellstring=movieindex{i};
        fprintf(fileID,'%s frames ',wellstring);
        frames=find(emptymatrix(i,:)>0);
        fprintf(fileID,'%d ',frames);
        fprintf(fileID,'\n');
    end
end
if sum(sum(blurredmatrix))>0
    fprintf(fileID,'blurred frames detected:\n');
    wellidx=find(sum(blurredmatrix,2)>0);
    for i=wellidx'
        wellstring=movieindex{i};
        fprintf(fileID,'%s frames ',wellstring);
        frames=find(blurredmatrix(i,:)>0);
        fprintf(fileID,'%d ',frames);
        fprintf(fileID,'\n');
    end
end
fclose(fileID);
type([datadir,'badframes.txt']);
save([datadir,'badframes.mat'],'movieindex','badframes');

cd([codepath,'Processing\']); %return to this directory