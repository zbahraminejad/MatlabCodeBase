function [traces1,traces2,mitoses] = combinedata(rowmat,colmat,sitemat,path,drugspike)
datadir = [path,'Data\'];
gatedir = [path,'GoodCells\'];
lastcommonmitosis = 1;      %default=1 (anything born after first post-drug mitosis might be a duplicate)
traces1 = [];
traces2 = [];
mitoses = [];
duplicates = 0;
for row=rowmat
    for col=colmat
        for site=sitemat
            moviename = [num2str(row),'_',num2str(col),'_',num2str(site)];
            load([datadir,moviename,'_alldata'],'best_rc');
            load([gatedir,moviename,'_goodcells'],'tracelist');
            load([gatedir,moviename,'_singletracedata'],'signal1','signal2','goodmitoses','realbirths');
            for cc=1:length(tracelist)
                i=tracelist(cc);
                eachmitosis = goodmitoses(i,:);
                postspikemitoses = sort(eachmitosis(eachmitosis>drugspike));
                if length(postspikemitoses)>lastcommonmitosis
                    postspikemitoses(1:lastcommonmitosis)=[];
                    if ismember(best_rc(i,1),postspikemitoses) && ismember(best_rc(i,2),tracelist) %'duplicate' data  
                        duplicates = duplicates+1;
                        tracelist(cc)=0;
                    end
                end
            end
            tracelist = tracelist(tracelist>0);
            traces1 = [traces1;signal1(tracelist,:)];
            traces2 = [traces2;signal2(tracelist,:)];
            mitoses = [mitoses;goodmitoses(tracelist,:)];
        end
    end
end
%fprintf(['duplicate count: ',num2str(duplicates),'\n']);