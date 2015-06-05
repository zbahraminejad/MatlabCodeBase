function sp=relocateD(data)
fmin=1;fmax=size(data,3);
sp=cell(size(data));

sp{fmin}=data{fmin};
x1=data{fmin}(:,1);y1=data{fmin}(:,2);
for f=fmin+1:fmax
    %fprintf('frame = %0.0f\n',f);
    sp{f}=zeros(size(sp{f-1}));
    %% preparation
    slo=size(x1,1);nlo=(x1>0);
    ser=(1:slo);ser=ser(nlo);
    x1=x1(nlo);x2=data{f}(:,1);X=double([x1;x2]);
    y1=y1(nlo);y2=data{f}(:,2);Y=double([y1;y2]);
    %z1=zeros(size(x1,1),1);z2=ones(size(x2,1),1);Z=[z1;z2];

    status12=zeros(size(x1,1),1);status21=zeros(size(x2,1),1);
    MMM=max(x1)-min(x1);
    dist12=ones(size(x1,1),1)*MMM;dist21=ones(size(x2,1),1)*MMM;
    neig12=cell(size(x1,1),1);neig12d=cell(size(x1,1),1);

    %% identifying static cells
    [NEWXY,ind]=sortrows([X,Y]);
    NEWXY_G=sum((NEWXY(2:end,:)-NEWXY(1:(end-1),:)).^2,2);
    NEWXY_lo=(NEWXY_G==0);NEWXY_lo1=logical([NEWXY_lo;0]);NEWXY_lo2=logical([0;NEWXY_lo]);
    stat1=ind(NEWXY_lo1);stat2=ind(NEWXY_lo2); % *sortrow will maintain initial order for equal values
    status12(stat1)=stat2-size(x1,1);dist12(stat1)=0;
    status21(stat2-size(x1,1))=stat1;dist21(stat2-size(x1,1))=0;

    oldlabel=(1:size(X,1))';statlabel=[stat1;stat2];
    X(statlabel)=[];Y(statlabel)=[];oldlabel(statlabel)=[];

    %% delaunay: pair closest cells from consecutive frames
    td=delaunay(X,Y);td_oldlabel=oldlabel(td);
    tdx=X(td);tdy=Y(td);
    tdl=single(td_oldlabel>size(x1,1));     %1 wherever cell from new frame

    if size(tdx,2)==1  %add as a fix for having few cells
        tdx=tdx';tdy=tdy';
    end

    X1=[tdx(:,1),tdx(:,1),tdx(:,2)];Y1=[tdy(:,1),tdy(:,1),tdy(:,2)];LO1=[tdl(:,1),tdl(:,1),tdl(:,2)];
    X2=[tdx(:,2),tdx(:,3),tdx(:,3)];Y2=[tdy(:,2),tdy(:,3),tdy(:,3)];LO2=[tdl(:,2),tdl(:,3),tdl(:,3)];
    DIS=sqrt((X2-X1).^2+(Y2-Y1).^2);LO=abs(LO2-LO1);
    DISA=sort(DIS,2);dis=median(DISA(:,1));
    %startp=[1,1,1,2,2,3];endp=[2,3,4,3,4,4];
    startp=[1,1,2];endp=[2,3,3];
    for i=1:size(td,1)
        temp_pos=td_oldlabel(i,:);
        temp_dis=DIS(i,:);
        temp_log=LO(i,:);
        for j=1:3
            if temp_log(j)>0        %cells are from different frames
                tempp=temp_pos([startp(j),endp(j)]);
                tp1=min(tempp);tp2=max(tempp)-size(x1,1);
                neig12{tp1}=[neig12{tp1},tp2]; %build list of neighbors from next frame
                neig12d{tp1}=[neig12d{tp1},temp_dis(j)]; %build list of distances from previous frame
                if temp_dis(j)<dist12(tp1)
                    dist12(tp1)=temp_dis(j);status12(tp1)=tp2;  %update nearest neighbor & distance from next frame
                end
                if temp_dis(j)<dist21(tp2)
                    dist21(tp2)=temp_dis(j);status21(tp2)=tp1;  %update nearest neighbor & distance from prev frame
                end
            end
        end
    end
    %% pair mutually assigned cells
    pair12=zeros(size(x1,1),1);pair21=zeros(size(x2,1),1);
    %nx2=zeros(size(x1,1),1);ny2=zeros(size(x1,1),1);
    for i=1:size(x2,1)
        if (status21(i)>0) && (status12(status21(i))==i) %cells from both frames mutually paired
            sp{f}(ser(status21(i)),:)=data{f}(i,:);      %save wellsss info from next frame to prev frame's cell id
            %nx2(status21(i))=x2(i);ny2(status21(i))=y2(i);
            pair21(i)=1;pair12(status21(i))=1;
        end
    end
    %% check remaining
    for i=(find(pair12==0))'                %not all unpaired cells from old frame are checked
        [temp_nei,m00]=unique(neig12{i});   %gather neighbors from new frame
        temp_dis=neig12d{i}(m00);           %gather associated distances
        temp_log=(pair21(temp_nei)==0)';    %are any of these neighbors from the new frame unpaired?
        if sum(temp_log)>0
            temp_nei=temp_nei(temp_log);    %gather the neighbors that are unpaired
            temp_dis=temp_dis(temp_log);    %gather assoc distances
            [B,IX]=min(temp_dis);           %find neighbor w/ minimum distance
            if B<dis            %less than median distance of all delaunay edges?
                %nx2(i)=x2(temp_nei(IX));ny2(i)=y2(temp_nei(IX));
                sp{f}(ser(i),:)=data{f}(temp_nei(IX),:);
                pair12(i)=1;pair21(temp_nei(IX))=1;
            end
        end
    end
    %nx2=[nx2;x2(~pair21)];ny2=[ny2;y2(~pair21)];
    sp{f}=[sp{f};data{f}(~pair21,:)];

    %% update values for next frame
    x1=sp{f}(:,1);y1=sp{f}(:,2);
end

end