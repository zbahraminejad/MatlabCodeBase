function [traces1,traces2,mitoses,firstunique] = combinedata_realbirths_history(rowmat,colmat,sitemat,path)
datadir = [path,'Data\'];
gatedir = [path,'GoodCells\'];
traces1 = [];
traces2 = [];
mitoses = [];
firstunique = [];
for row=rowmat
    for col=colmat
        for site=sitemat
            %moviename = [num2str(row),'_',num2str(col),'_',num2str(site)];
            moviename = wellnum2str(row,col,site);
            load([datadir,moviename,'_alldata_IF'],'best_rc');
            load([gatedir,moviename,'_goodcells'],'tracelist');
            load([gatedir,moviename,'_singletracedata'],'signal1','signal2','goodmitoses','realbirths');
            for cc=1:length(tracelist)
                i=tracelist(cc);
                eachmitosis = goodmitoses(i,:);
                eachmitosis = sort(unique(eachmitosis(eachmitosis>=realbirths(i))));
                goodmitoses(i,:) = 0;
                goodmitoses(i,1:length(eachmitosis)) = eachmitosis;
            end
            traces1 = [traces1;signal1(tracelist,:)];
            traces2 = [traces2;signal2(tracelist,:)];
            mitoses = [mitoses;goodmitoses(tracelist,:)];
            firstunique = [firstunique;realbirths(tracelist,:)];
        end
    end
end