function Ancestors = getAncestors(ChildInd, Ancestors, EdgeList)
    ParentIndList = find(~cellfun(@isempty,EdgeList(:,ChildInd)))';
    for i=1:size(ParentIndList,2)
        Ancestors = union(Ancestors,EdgeList{ParentIndList(i),ChildInd});
        newChildInd = ParentIndList(i);
        %newParentIndList = find(~cellfun(@isempty,EdgeList{:,newChildInd}))';
        Ancestors = getAncestors(newChildInd, Ancestors, EdgeList);
    end
end