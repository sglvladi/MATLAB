DataList = {};
timestamps = [];
i=1;
j=1;
while i<=size(radar,1)
    tempDataList = [];
    timestamp = radar(i,7);
    z = 1;
    if isnan(timestamp)
        i = i+1;
        continue;
    end
   
    while(i<=size(radar,1)&&timestamp==radar(i,7))
       tempDataList(z,:) = radar(i,3:4)/1000;
       i=i+1
       z=z+1;
    end
    DataList{j} = tempDataList';
    timestamps(j) =timestamp;
    j=j+1
end

Dt = timestamps(2:end) - timestamps(1:end-1);
i=1;
while i<=size(Dt,2)
    if(Dt(i)==0)
        DataList{i-1} = [DataList{i-1},DataList{i}];
        DataList(i) = [];
        Dt(i)=[];
    else
        i = i+1;
    end
end


% for i=1:size(DataList,2)
% plot(DataList{i}(:,1), DataList{i}(:,2), 'r*');
% axis([-2000 -200 -3000 3000]);
% pause(0.1);
% end