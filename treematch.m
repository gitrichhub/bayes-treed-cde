function out = treematch(T1,T2,nodeId1,nodeId2)
    out = 1;
    % Find root node if empty nodes Ids specified
    if isempty(nodeId1)
        for jj = 1:length(T1.Allnodes)
            if isempty(T1.Allnodes{jj}.Parent)
                nodeId1 = T1.Allnodes{jj}.Id;
            end
        end
    end
    if isempty(nodeId2)
        for jj = 1:length(T2.Allnodes)
            if isempty(T2.Allnodes{jj}.Parent)
                nodeId2 = T2.Allnodes{jj}.Id;
            end
        end
    end


    nind1 = nodeind(T1,nodeId1);
    nind2 = nodeind(T2,nodeId2);
    node1 = T1.Allnodes{nind1};
    node2 = T2.Allnodes{nind2};
    if ~isempty(node1.Rule) && ~isempty(node2.Rule)
        if node1.Rule{1} == node2.Rule{1}
            if isa(node1.Rule{2},'cell')
                rule1 = sort(node1.Rule{2});
                rule2 = sort(node2.Rule{2});
                if ~all(strcmp(rule1,rule2))
                    out = 0;
                    return;
                end
            else
                if node1.Rule{2} ~= node2.Rule{2}
                    out = 0;
                    return;
                end
            end
            % If the rules match, recursively check the Left and right children
            % children
            lchild1 = node1.Lchild;
            lchild2 = node2.Lchild;
            nindlchild1 = nodeind(T1,lchild1);
            nindlchild2 = nodeind(T2,lchild2);
            nodelchild1 = T1.Allnodes{nindlchild1};
            nodelchild2 = T2.Allnodes{nindlchild2};
            if isempty(nodelchild1.Rule) && isempty(nodelchild2.Rule) % both empty: ok
             elseif isempty(nodelchild1.Rule)  || isempty(nodelchild2.Rule) % one empty: don't match
                out = 0;
                return;
            else % both not empty
                out = treematch(T1,T2,lchild1,lchild2);
            end

            if out
                rchild1 = node1.Rchild;
                rchild2 = node2.Rchild;
                nindrchild1 = nodeind(T1,rchild1);
                nindrchild2 = nodeind(T2,rchild2);
                noderchild1 = T1.Allnodes{nindrchild1};
                noderchild2 = T2.Allnodes{nindrchild2};
                if isempty(noderchild1.Rule) && isempty(noderchild2.Rule) % both empty: ok
                elseif isempty(noderchild1.Rule)  || isempty(noderchild2.Rule) % one empty: don't match
                    out = 0;
                    return;
                else % both not empty
                    out = treematch(T1,T2,rchild1,rchild2);
                end
                if ~out
                    return;
                end
            else
                return;
            end
        else
            out = 0; 
            return;
        end
    elseif any([~isempty(node1.Rule), ~isempty(node2.Rule)]) % they don't match
        out = 0;
        return;
    end
end