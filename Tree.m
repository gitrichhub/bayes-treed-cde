classdef Tree
    properties
        % X
        Allnodes % array of all Nodes
        NodeIds % A vector with the Ids of the nodes in Allnodes
        Xclass % column types
        ncol
        Emptynodes
        Varnames
        GP % Default GP structure
        ygrid
        ygridn
        m % grid points for y
        Lliketree % Tree log-likelihood
        %Leafmin % minimum number of values at a leaf
    end
    methods
        % Constructor
        function out = Tree(y,X)
            % Store data
            % out.X = X;
            % Create root node
            rootnode = Nodes(0,[],[],[],[],1:size(X,1));
            out.Allnodes{1} = rootnode;
%             if ~isempty(Leafmin)
%                 out.Leafmin = Leafmin;
%             end
            out.NodeIds = 0;
            out.Emptynodes = 0;
            out.Varnames = X.Properties.VariableNames;
            % ygrid
            out.m = 400;
            m = out.m;
            ymin = min([min(y), mean(y) - 3*std(y)]);
            ymax = max([max(y), mean(y) + 3*std(y)]);
            ygrid = linspace(ymin,ymax,m);
            ygridn = zscore(ygrid); % normalized ygrid;
            out.ygrid = ygrid;
            out.ygridn = ygridn;
            
            
            % Get column types
            if isa(X,'table')
                out.ncol = size(X,2);
                vartypes = varfun(@class,X,'OutputFormat','cell');
                out.Xclass = vartypes;
                if ~all(strcmp('double',vartypes) | strcmp('cell',vartypes))
                    error('Data columns of X must be "double" or "cell".')
                end
            else
                error('Design matrix must be of class "table"');
            end
            % GP structure
            % Set up the general GP structure
            % Priors for the sigma2 and l (Assuming the grid points have been centered
            %   and scaled)
            pm = prior_sqrtt('s2',10,'nu',4);
            pl = prior_t('s2', 1, 'nu', 4);
            % Defautl smoothing values
            sigma2 = 1;
            % Best guess of lengthscale (From gpsmooth function in lgpdens from
            %    riihimaki's excellent MATLAB code)
            %Xn = zscore(X);
            h=max(diff(ygridn(end,1:2)).^2,1/length(y).^(1/5)/2);
            % With Prior
            % NOTE: the starting value for lengthScale was obtained by looking in the
            %   gpsmooth function in the lgpdens function from Riihimaki.
            cf = gpcf_sexp('magnSigma2',sigma2,'magnSigma2_prior',pm,...
                'lengthScale_prior',pl,'lengthScale',h*repmat(2,[1 size(ygridn,1)]));
            gpmflin = gpmf_linear('prior_mean',0,'prior_cov',100);
            %gpmfsq = gpmf_squared('prior_mean',0,'prior_cov',100,'interactions','on');
            % NOTE: I have turned interactions off since this was not mentioned in
            %    the original LGP density paper for 1D (probably doesn't matter for 1D).
            gpmfsq = gpmf_squared('prior_mean',0,'prior_cov',100,'interactions','off');
            % Set GP
            out.GP = gp_set('lik',lik_lgp,'cf',cf,'jitterSigma2',1e-6,'meanf',{gpmflin,gpmfsq},...
                'latent_method', 'Laplace');
            
            % Calculate log-likelihood of root node
            out = llike(out,0,y);
            % Set tree log-likelihood to that of the root node;
            out.Lliketree = out.Allnodes{1}.Llike;
        end
        
        % Calculate log-likelihood of a node
        function out = llike(obj,nodename,y)
            nind = nodeind(obj,nodename);
            out = obj;
            thenode = obj.Allnodes{nind};
            thenode = loglikefunc(thenode,obj,y);
            thenode.Updatellike = 0; 
            out.Allnodes{nind} = thenode;
        end
            
        % Birth Function
        function out = birth(obj,y,X)
            out = obj;
            % Get Terminal Node indices and IDs
            [Ind,tnodeIDs] = termnodes(obj);
            % Randomly select a terminal node
            RIND = randsample(length(Ind),length(Ind));
            
            for ii = 1:length(RIND)
                rind = RIND(ii);
                birthID = tnodeIDs(rind);
                birthindex = Ind(rind);

                % Birth at the terminal node;
                % Find indices of data at current node;
                % Xind = nodeData(obj,birthindex);
                % Determine Split Rule based on data

                node = obj.Allnodes{birthindex};
                nrule = drawrule(node,obj,X,0);
                if ~isempty(nrule)
                    % Assign new rule
                    node = newrule(node,nrule);
                    out.Allnodes{birthindex} = node;  

                    % Add children nodes;
                    % Find unique new IDs for children
                    newids = newIDs(obj);
                    % Get data which is being passed down from parent
                    [XindL,XindR] = childrendata(out,birthID,X); % must use out since it has the new rule
                    % Add Nodes
                    childnode1 = Nodes(newids(1),birthID,[],[],[],XindL);
                    childnode2 = Nodes(newids(2),birthID,[],[],[],XindR);
                    out = addnode(out,childnode1,birthindex,'L',y);
                    out = addnode(out,childnode2,birthindex,'R',y);
                    % Update log-likelihood
                    nn = nnodes(out);
                    out.Lliketree = out.Lliketree - node.Llike + ...
                        out.Allnodes{nn-1}.Llike + out.Allnodes{nn}.Llike;
                    
                    return;
                end
            end
            error('Birth Step Not Possible')
        end
        
        % Prune Function
        function out = prune(obj,y)
            if length(obj.Allnodes) > 1
                % Find parents of two-terminal nodes
                [I,~] = terminalparents(obj);
                % Uniformly choose one node to be pruned
                if length(I) > 1
                    pind = randsample(I,1);
                else
                    pind = I;
                end
                % Prune the node
                out = delnode(obj,pind,y);              
            else
                error('Cannot prune a root node.')
            end
        end
        
        % Change Function
        function out = change(obj,y,X)
            % Find interior nodes
            [I,~] = interiornodes(obj);
            if isempty(I)
                error('Cannot perform a change step on a root node with no children.')
            end
            
            % Generate a new split rule from the prior (CART)
            notok = 1;
            cntr = 0;
            maxiter = 100;
            colind = 0; % Choose a random column in drawrule function for first time;
            out = obj;
            while notok
                % Randomly select an interior node to change the rule
                changeind = randsample(I,1);
                
                % Reset on each iteration...
                node = obj.Allnodes{changeind};
                out = obj;
                [nrule, colind] = drawrule(node,obj,X,colind);
                node = newrule(node,nrule);
                out.Allnodes{changeind} = node;
                out.Emptynodes = 0;
                out = descendentdata(out,node.Id,X);
                if out.Emptynodes == 0
                    notok = 0;
                end
                cntr = cntr + 1;
                if cntr > maxiter
                    error('Change step could not find a satisfactory rule in 100 iterations')
                end               
            end
            % Update log-likelihood on terminal nodes which need updating
            % Also update the Tree log-likelihood
            out.Lliketree = 0;
            [I,Ids] = termnodes(out);
            for ii=1:length(Ids)
                if out.Allnodes{I(ii)}.Updatellike
                    out = llike(out,Ids(ii),y);
                end
                out.Lliketree = out.Lliketree + out.Allnodes{I(ii)}.Llike;
            end
        end
        
        % Swap
        function out = swap(obj,y,X)
            % Find internal nodes whose parent is also an internal node;
            Ids = parentchildpairs(obj);
            if isempty(Ids)
                error('Cannot perform swap step: No internal parent-child pairs.')
            end
            Ids = randsample(Ids,length(Ids));
            % Loop through until we find a swap without empty nodes.
            for ii = 1:length(Ids)
                tmpout = obj;
                % Randomly select one child
                rind = Ids(ii);
                nind = nodeind(obj,rind);
                node = obj.Allnodes{nind};
                % Find the node's parent
                nindParent = nodeind(obj,node.Parent);
                nodeParent = obj.Allnodes{nindParent};
                % Parent Rule
                prule = nodeParent.Rule;
                % Is the node a left or right child?
                LRind = ismember([nodeParent.Lchild,nodeParent.Rchild],rind);
                indR = nodeind(obj,nodeParent.Rchild);
                indL = nodeind(obj,nodeParent.Lchild);
                if all(LRind == [1 0]) % Left Child
                    lrule = node.Rule;
                    %indR = nodeind(obj,nodeParent.Rchild);
                    rrule = obj.Allnodes{indR}.Rule;
                elseif all(LRind == [0 1]) % Right Child
                    rrule = node.Rule;
                    %indL = nodeind(obj,nodeParent.Lchild);
                    lrule = obj.Allnodes{indL}.Rule;
                end

                % Are the two children rules equal?
                eqrule = SplitRule.equalrule(lrule,rrule);
                if eqrule
                    % Switch Rules
                    newprule = lrule;
                    newLrule = prule;
                    newRrule = prule;
                else
                    if all(LRind == [1 0]) % Left Child
                        newprule = lrule;
                        newLrule = prule;
                        newRrule = rrule;
                    elseif all(LRind == [0 1]) % Right Child
                        newprule = rrule;
                        newLrule = lrule;
                        newRrule = prule;
                    end
                end
                tmpout.Allnodes{nindParent}.Rule = newprule;
                tmpout.Allnodes{indL}.Rule = newLrule;
                tmpout.Allnodes{indR}.Rule = newRrule;

                % Update Descendent Data
                tmpout.Emptynodes = 0;
                tmpout = descendentdata(tmpout,nodeParent.Id,X);
                if tmpout.Emptynodes == 0
                    out = tmpout;
                    % Update log-likelihood on terminal nodes which need updating
                    % Also update Lliketree
                    out.Lliketree = 0;
                    [I,Ids] = termnodes(out);
                    for jj=1:length(Ids)
                        if out.Allnodes{I(jj)}.Updatellike
                             out = llike(out,Ids(jj),y);
                        end
                        out.Lliketree = out.Lliketree + out.Allnodes{I(jj)}.Llike;
                    end
                    return;
                end
            end
            error('Swap step not possible.')
            
        end
            
              
        % Find interior nodes (Ids) with a parent who is also an interior node;
        function Ids = parentchildpairs(obj)
            [~,Id] = termnodes(obj);
            nIds = obj.NodeIds;
            nIds = nIds(nIds > 0); % Exclude the root node;
            ind = ~ismember(nIds,Id); 
            Ids = nIds(ind);
            if min(size(Ids)) == 0
                Ids = [];
            end
        end
       
        % delete node (get rid of children and split rule)
        function out = delnode(obj,nodeindex,y)
           out = obj;
           node = obj.Allnodes{nodeindex};
           % Delete children (Somebody help them!)
           indL = nodeind(obj,node.Lchild);
           indR = nodeind(obj,node.Rchild);
           
           % Subtract log-likelihood contributions of deleted children
           out.Lliketree = out.Lliketree - out.Allnodes{indL}.Llike -...
               out.Allnodes{indR}.Llike;           
           
           out.Allnodes([indL,indR]) = [];
           %out.Allnodes(indR) = [];
           out.NodeIds([indL, indR]) = [];
           %out.NodeIds(indR) = [];
           % Get rid of children on parent node
           node.Lchild = [];
           node.Rchild = [];
           node = newrule(node,[]);
           out.Allnodes{nodeindex} = node; % update node
           % Update likelihood if necessary
           if node.Updatellike
               out = llike(out,out.Allnodes{nodeindex}.Id,y);
           end
           % Add in Log-likelihood of the parent node
           out.Lliketree = out.Lliketree + out.Allnodes{nodeindex}.Llike;
        end
        
        % Add node to tree
        function out = addnode(obj,node,parentind,LR,y)
            nn = nnodes(obj);
            out = obj;
            out.Allnodes{nn + 1} = node;
            out.NodeIds(nn + 1) = node.Id;
            if strcmp(LR,'L')
                out.Allnodes{parentind}.Lchild = node.Id;
            elseif strcmp(LR,'R')
                out.Allnodes{parentind}.Rchild = node.Id;
            else
                error('LR must be either "L" or "R"');
            end
            % Calculate log-likelihood of new node
            out = llike(out,node.Id,y);
        end
        
        % Returns the index and Ids of the interior nodes
        function [I,Id] = interiornodes(obj)
            nn = nnodes(obj);
            if nn <= 1
               I = [];
               Id = [];
               return;
            end
            Id = zeros(nn,1);
            I = Id;
            cntr = 1;
            for ii=1:nn
                if ~isempty(obj.Allnodes{ii}.Rule)
                    Id(cntr) = obj.NodeIds(ii);
                    I(cntr) = ii;
                    cntr = cntr + 1;
                end
            end
            Id = Id(1:(cntr - 1));
            I = I(1:(cntr - 1));          
        end

        
        % Returns the index and Ids of the terminal nodes
        function [I,Id] = termnodes(obj)
            nn = nnodes(obj);
            Id = zeros(nn,1);
            I = Id;
            cntr = 1;
            for ii=1:nn
                if isempty(obj.Allnodes{ii}.Rule)
                    Id(cntr) = obj.NodeIds(ii);
                    I(cntr) = ii;
                    cntr = cntr + 1;
                end
            end
            Id = Id(1:(cntr - 1));
            I = I(1:(cntr - 1));
        end
        
        % Gets new unique IDs for birth step
        function out = newIDs(obj)
            nn = nnodes(obj);
            %allIds = zeros(nn,1);
            %for ii=1:nn
            %    allIds(nn) = obj.Allnodes{ii}.Id;
            %end
            allIds = obj.NodeIds;
            alln = 0:(nn-1);
            I = ~ismember(allIds,alln);
            out = alln(I);
            cntr = 0;
            while length(out) < 2
                out = [out,nn + cntr];
                cntr = cntr + 1;
            end
        end
        
        % Returns the number of nodes (including leaves and root)
        function out = nnodes(obj)
           out = length(obj.Allnodes); 
        end
        
        % Obtain the data for each of the children of a designated node
        function [XindL,XindR,empty] = childrendata(obj,nodeid,X)
            nind = nodeind(obj,nodeid);
            node = obj.Allnodes{nind};
            therule = node.Rule; % SplitRule class object
            if isempty(therule)
                error('No rule specified for the node.')
            end
            if strcmp(therule.Varclass,'double')
                %I = find(table2array(X(node.Xind,therule.Varcol)) <= therule.Varrule);
                I = table2array(X(node.Xind,therule.Varcol)) <= therule.Varrule ;
            elseif strcmp(therule.Varclass,'cell')
                I = ismember(table2cell(X(node.Xind,therule.Varcol)),therule.Varrule);
            else
                error('Unexpected variable class found.')
            end
            XindL = node.Xind(I);
            XindR = node.Xind(~I); 
            
            % Check for empty children (children with no data)
            emptyL = 0;
            emptyR = 0;
            if min(size(XindL)) == 0
                emptyL = 1;
                XindL = [];
            end
            if min(size(XindR)) == 0
                emptyR = 1;
                XindR = [];
            end
            empty = [emptyL,emptyR];               
            
        end
        
        % Update data of all descendents of a node (recursive function)
        function out = descendentdata(obj,nodeid,X)
            nind = nodeind(obj,nodeid);
            node2 = obj.Allnodes{nind};
            if ~isempty(node2.Lchild) && ~isempty(node2.Rchild)
                [XindL,XindR,empty] = childrendata(obj,nodeid,X);           
                out = obj;
                if any(empty)
                    out.Emptynodes = 1;
                end
                nindL = nodeind(out,node2.Lchild);
                nindR = nodeind(out,node2.Rchild);
                out = updatedata(out,nindL,XindL);
                out = updatedata(out,nindR,XindR);
                % Now do it for the descendents
                out = descendentdata(out,node2.Lchild,X);
                out = descendentdata(out,node2.Rchild,X);
            else
                out = obj;
            end
        end
        
        
        
        % Get the index of a node with nodeid;
        function out = nodeind(obj,nodeid)
            nn = length(nodeid);
            if nn == 1
                out = find(obj.NodeIds == nodeid);
            elseif nn > 1
                out = zeros(nn,1);
                for ii=1:nn
                    out(ii) = find(obj.NodeIds == nodeid(ii));
                end
            else
                error('Something went awry')
            end
        end
        
        % Find parents of children who are BOTH terminal nodes;
        function [I,ids] = terminalparents(obj)
            % Find Term Nodes
            [I,~] = termnodes(obj);
            % If only the root node, return the root node.
            if length(I) < 2 && length(obj.Allnodes) < 2
                I = 1;
                ids = obj.Allnodes{I}.Id;
                return
            end

            % Find parents of terminal nodes
            parentIDs = zeros(length(I),1);
            for ii = 1:length(I)
                parentIDs(ii) = obj.Allnodes{I(ii)}.Parent;
            end
            if length(parentIDs) < 1
               error('There are no nodes with two terminal children.') 
            end
            % Find which parentIDs have two terminal node children
            [C,~,IC] = unique(parentIDs);
            tab = tabulate(C(IC));
            tabind = tab(:,2) > 1;
            ids = tab(tabind,1);
            I = nodeind(obj,ids);
        end
        
        % Update the data index of the node at index nodeind
        function out = updatedata(obj,nodeind,Xind)
            out = obj;
            node1 = out.Allnodes{nodeind};
            % Determine if the data has changed and mark it.
            if length(node1.Xind) ~= length(Xind) % Different sizes
                node1.Updatellike = 1;
            elseif ~all(sort(Xind) == sort(node1.Xind))
                node1.Updatellike = 1;
            end
            
            node1.Xind = Xind;
            out.Allnodes{nodeind} = node1;
        end
            
        
        % Returns the index of X which reach the node specified
        function out = nodeData(obj,nodeindex)
%             for ii = 1:length(nnodes(obj))
%                 if obj.Allnodes{ii}.Id == nodeID
%                     node = obj.Allnodes{ii};
%                     break
%                 end
%             end
            node = obj.Allnodes{nodeindex};
            if exist('node','var')
                if isempty(node.Lchild) && isempty(node.Rchild)
                    out = node.Xind;
                else
                    error('You have not yet written code for this...')
                end 
            else
                error('No node with that ID found')
            end
        end
        
        % Graphs
        function treelines(obj,nodename,level,treedepth,parentxloc,LR)
            width = 1; % space between terminal nodes
            maxloc = 1 + width*(2*treedepth - 1);
            nind = nodeind(obj,nodename);
            node = obj.Allnodes{nind};
            if nodename > 0
                if LR == 'L'
                    plusminus = -1;
                elseif LR == 'R'
                    plusminus = 1;
                else
                    error('Must correctly specify "LR" as "L" or "R"')
                end
                delta = width*2^(treedepth + level - 1);
                %delta = width - 1/(2*treedepth);
                %delta = (1 + width*(2*treedepth - 1) - (maxloc+1)/2)/treedepth;
                xval = parentxloc + plusminus*delta;
                plot([parentxloc, xval],[level,level-1],'b')
                
                
                               
            else
                % Get the xval to pass to children
                %maxloc = 2*treedepth*width;     
                xval = (maxloc + 1)/2;
                %parentxloc = xval; % just for the root node...
            end
            
            % Plot the rule of the parent
            if ~isempty(node.Rule) % if parent node
                
                colnum = node.Rule.Varcol;
                colname = obj.Varnames(colnum);
                if strcmp(node.Rule.Varclass,'double')
                    ruletext = strcat(colname,' \leq ',num2str(node.Rule.Varrule));
                else
                    for ii=1:length(node.Rule.Varrule)
                        if ii == 1
                            grp = node.Rule.Varrule{ii};
                        else
                            grp = strcat(grp,' , ',node.Rule.Varrule{ii});
                        end
                    end
                    ruletext = strcat(colname,' \in \{',grp,'\}');
                    %ruletext = strcat(colname,' \in \{',cell2mat(node.Rule.Varrule),'\}');
                end
                text(xval,level-1,ruletext,'HorizontalAlignment','right')
            end
            
            %nn = length(obj.Allnodes);
            if ~isempty(node.Lchild) && ~isempty(node.Rchild)
                treelines(obj,node.Lchild,level-1,treedepth,xval,'L')
                treelines(obj,node.Rchild,level-1,treedepth,xval,'R')
            end
        end
        
        
        function Treeplot(obj)
            nn = length(obj.Allnodes);
            if nn <= 1
                error('Tree must be more than a root node.')
            else
                nodegraph = zeros(1,nn);
                for ii = 2:nn
                   nodegraph(ii) = obj.Allnodes{ii}.Parent + 1;
                end
                %treeplot(nodegraph)
                treedepth = max(nodegraph);
                %width = 1; % space between terminal nodes
                %maxloc = 2*treedepth;
                figure;
                hold on;
                treelines(obj,0,0,treedepth,'','')
                hold off;
                axis off;
            end
        end
    end
end