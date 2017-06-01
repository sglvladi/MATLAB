function D = drawTreeObj(TreeObj, type, plot_title)
    if strcmp('MO',type)
        figure;
        % Draw Tree 
        % Input: TreeObj
        source = [];
        target = [];
        names = {};
        remainders = {};
        weights = {};
        for NodeInd=1:size(TreeObj.NodeList,2)
            % Generate source and target arrays
            NodeObj = TreeObj.NodeList{1, NodeInd};
            for i = 1:size(NodeObj.ChildIndList,2)
                ChildInd = NodeObj.ChildIndList(i);
                source = [source, NodeInd];
                target = [target, ChildInd];
                %Generate edge weights
                weight = strcat(weight, sprintf('%d,', TreeObj.EdgeList(NodeInd, ChildInd)));
                weights{end+1} = weight;
            end
            % Generate node names
            if(NodeInd==1)
                name = '$$(\emptyset, (';  
            else
                name = sprintf('$$(M%d, (', NodeObj.MeasInd);
            end
            if ~size(NodeObj.TrackIndList,2)
               name = strcat(name, '\emptyset'); 
            else
                for i= 1:size(NodeObj.TrackIndList,2)
                    TrackInd = NodeObj.TrackIndList(i);
                    if i==size(NodeObj.TrackIndList,2)
                        name = sprintf(strcat(name, sprintf('T%d',TrackInd)));
                    else
                        name = sprintf(strcat(name, sprintf('T%d,',TrackInd)));
                    end
                        %name = sprintf('%sT%d,', [name, TrackInd]);
                end 
            end
            name = strcat(name, sprintf(') )$$'));
            remainders_tmp = '$$[';
            if ~size(NodeObj.Remainders,1)
               remainders_tmp = strcat(remainders_tmp, '\emptyset'); 
            else
                for i = 1:size(NodeObj.Remainders,1)
                   TrackInd = NodeObj.Remainders(i,1);
                    if i==size(NodeObj.Remainders,1)
                        remainders_tmp = strcat(remainders_tmp, sprintf('T%d', TrackInd));
                    else
                        remainders_tmp = strcat(remainders_tmp, sprintf('T%d,', TrackInd));
                    end
                end
            end
            remainders_tmp = strcat(remainders_tmp, ']$$');
            names{NodeInd} = name;
            remainders{NodeInd} = remainders_tmp;
        end
        celldata = cellstr(names);
        D = digraph(source,target);%, randi(100, size(source)), celldata);
        D.Edges.Power = weights';

        P = plot(D, 'Layout', 'layered', 'EdgeLabel',D.Edges.Power);

        for i=1:size(names,2)
             txt = sprintf('%s \newline %s ',  [names{i}, remainders{i}]);
             name = names{i};
             remainder = remainders{i};
             text(P.XData(i)+0.05,P.YData(i)+0.05-0.05,{name,remainder},'HorizontalAlignment','left', 'interpreter','latex');
        end
    else
        figure;
        % Draw Net 
        % Input: TreeObj
        source = [];
        target = [];
        names = {};
        remainders = {};
        weights = {};
        for NodeInd=1:size(TreeObj.NodeList,2)
            % Generate source and target arrays
            NodeObj = TreeObj.NodeList{1, NodeInd};
            for i = 1:size(NodeObj.ChildIndList,2)
                ChildInd = NodeObj.ChildIndList(i);
                source = [source, NodeInd];
                target = [target, ChildInd];
                %Generate edge weights
                weight = '';
                weight = strcat(weight, sprintf('%d', full(TreeObj.EdgeList(NodeInd, ChildInd))));   
                weights{end+1} = weight;
            end
            % Generate node names
            if(NodeInd==1)
                name = '$$(\emptyset, (';  
            else
                name = sprintf('$$(T%d, (', NodeObj.TrackInd);
            end
            if ~size(NodeObj.MeasInd,2)
               name = strcat(name, '\emptyset '); 
            else
               name = sprintf(strcat(name, sprintf('M%d',NodeObj.MeasInd)));
            end
            name = strcat(name, sprintf(') )$$'));
            remainders_tmp = '$$[';

            if ~size(NodeObj.Remainders,1)
               remainders_tmp = strcat(remainders_tmp, '\emptyset'); 
            else
                for i = 1:size(NodeObj.Remainders,1)
                   TrackInd = NodeObj.Remainders(i,1);
                    if i==size(NodeObj.Remainders,1)
                        remainders_tmp = strcat(remainders_tmp, sprintf('M%d', TrackInd));
                    else
                        remainders_tmp = strcat(remainders_tmp, sprintf('M%d,', TrackInd));
                    end
                end
            end
            remainders_tmp = strcat(remainders_tmp, ']$$');
            names{NodeInd} = name;
            remainders{NodeInd} = remainders_tmp;
        end
        celldata = cellstr(names);
        D = digraph(source,target);%, randi(100, size(source)), celldata);
        D.Edges.Power = weights';

        P = plot(D, 'Layout', 'layered', 'EdgeLabel',D.Edges.Power);
%         highlight(P,3,5,'EdgeColor','red')
%         highlight(P,4,6,'EdgeColor','red')
%         highlight(P,6,7,'EdgeColor','red')
        for i=1:size(names,2)
             txt = sprintf('%s \newline %s ',  [names{i}, remainders{i}]);
             name = names{i};
             remainder = remainders{i};
             text(P.XData(i)+0.05,P.YData(i)+0.05-0.05,{name},'HorizontalAlignment','left', 'interpreter','latex');
        end
    end
    title(plot_title)
end