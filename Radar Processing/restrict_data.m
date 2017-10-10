i=1;
bounding_box = [-2 -.800 2 3];% [-.700 -.400 -.700 .400];
while i<=size(DataList,2)
    j=1;
    while j<=size(DataList{i},2)
        if DataList{i}(1,j)>bounding_box(2) || DataList{i}(1,j)<bounding_box(1) || DataList{i}(2,j)>bounding_box(4) || DataList{i}(2,j)<bounding_box(3)
            DataList{i}(:,j) = [];
        else
            j=j+1;
        end
    end
    if size(DataList{i},2)==0 && i<size(DataList,2)-1
        DataList(i) = [];
        i
        Dt(i+1) = Dt(i+1) + Dt(i);
        Dt(i) = [];
    else
        i = i+1;
    end
end