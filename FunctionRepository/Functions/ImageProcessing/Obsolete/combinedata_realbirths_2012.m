function [traces1,traces2,mitoses] = combinedata_realbirths(rowmat,colmat,sitemat,path)
datadir = [path,'Data\'];
gatedir = [path,'GoodCells\'];
traces1 = [];
traces2 = [];
mitoses = [];
for row=rowmat
    for col=colmat
        for site=sitemat
            moviename = [num2str(row),'_',num2str(col),'_',num2str(site)];
            load([datadir,moviename,'_alldata_HS_L05'],'best_rc');
            load([gatedir,moviename,'_goodcells'],'tracelist');
            load([gatedir,moviename,'_singletracedata'],'signal1','signal2','goodmitoses','realbirths');
            for cc=1:length(tracelist)
                i=tracelist(cc);
                eachmitosis = goodmitoses(i,:);
                eachmitosis = sort(eachmitosis(eachmitosis>=realbirths(i)));
                goodmitoses(i,:) = 0;
                goodmitoses(i,1:length(eachmitosis)) = eachmitosis;
                if realbirths(i)>1
                    signal1(i,1:realbirths(i)-1) = -10000;  %ensures that only unique data is analyzed
                    signal2(i,1:realbirths(i)-1) = -10000;
                end
            end
            traces1 = [traces1;signal1(tracelist,:)];
            traces2 = [traces2;signal2(tracelist,:)];
            mitoses = [mitoses;goodmitoses(tracelist,:)];
        end
    end
end