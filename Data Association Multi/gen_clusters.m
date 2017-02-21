clusters = {}

for i=1:size(ValidationMatrix,1) % Iterate over all measurements
    matched =[];
    temp_clust = find(ValidationMatrix(i,:))'; % Extract associated tracks
    if (~isempty(temp_clust))   % If measurement matched with any tracks
        % Check if matched tracks are members of any clusters
        for j=1:size(clusters,2)
            a = ismember(temp_clust, cell2mat(clusters(1,j)));
            if (ismember(1,a)~=0)
                matched = [matched, j]; % Store matched cluster ids
            end   
        end
        if(size(matched,2)==1) % If only matched with a single cluster, join.
            clusters{1,matched(1)}=union(cell2mat(clusters(1,matched(1))), temp_clust);
        elseif (size(matched,2)>1) % If matched with more that one clusters
            matched = sort(matched); % Sort cluster ids
            % Start from last cluster, joining each one with the previous
            %   and removing the former.  
            for match_ind = size(matched,2)-1:-1:1
                clusters{1,match_ind}=union(cell2mat(clusters(1,match_ind)), cell2mat(clusters(1,match_ind+1)));
                clusters(:,match_ind+1)=[];
            end
            % Finally, join with associated track.
            clusters{1,match_ind}=union(cell2mat(clusters(1,match_ind)), temp_clust);
        else % If not matched with any cluster, then create a new one.
            clusters{1,size(clusters,2)+1} = temp_clust;
        end
    end
end