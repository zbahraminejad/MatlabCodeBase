function [traces1,traces2,mitoses,firstunique,EdUlist,DAPIlist] = combinedata_EdU(rowmat,colmat,sitemat,path)
datadir = [path,'Data\'];
gatedir = [path,'GoodCells\'];
traces1 = [];
traces2 = [];
mitoses = [];
firstunique = [];
EdUlist = [];
DAPIlist = [];
for row=rowmat
    for col=colmat
        for site=sitemat
            moviename = wellnum2str(row,col,site);
            %moviename = [num2str(row),'_',num2str(col),'_',num2str(site)];
            load([datadir,moviename,'_alldata'],'bestsp','best_rc');
            load([gatedir,moviename,'_goodcells'],'tracelist');
            load([gatedir,moviename,'_singletracedata'],'signal1','signal2','goodmitoses','realbirths');
            for cc=1:length(tracelist)
                i=tracelist(cc);
                eachmitosis = goodmitoses(i,:);
                eachmitosis = sort(unique(eachmitosis(eachmitosis>=realbirths(i))));
                goodmitoses(i,:) = 0;
                goodmitoses(i,1:length(eachmitosis)) = eachmitosis;
            end
            %%% add EdU and DAPI data %%%%%%%%%%
            totalframes = size(bestsp,3);
            EdUstain = bestsp{totalframes}(tracelist,8);
            DAPIstain = bestsp{totalframes}(tracelist,9).*bestsp{totalframes}(tracelist,4);
            %%% concatenate master list %%%%%%%%
            traces1 = [traces1;signal1(tracelist,:)];
            traces2 = [traces2;signal2(tracelist,:)];
            mitoses = [mitoses;goodmitoses(tracelist,:)];
            firstunique = [firstunique;realbirths(tracelist,:)];
            EdUlist = [EdUlist;EdUstain];
            DAPIlist = [DAPIlist;DAPIstain];
        end
    end
end