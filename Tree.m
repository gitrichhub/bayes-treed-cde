classdef Tree
    % TODO: remove parentchildagree functions...
    % TODO: Store tree prior as a class property.
    properties
        % X
        Allnodes % array of all Nodes
        NodeIds % A vector with the Ids of the nodes in Allnodes
        Xclass % column types
        ncol
        Smallnodes
        Varnames
        GP % Default GP structure
        ygrid
        ygridn
        m % grid points for y
        Lliketree % Tree log-likelihood (not computed with temp)      
        Leafmin % minimum number of values at a leaf
        gamma
        beta
        Prior
        Temp
        Ntermnodes
        %Updateprior        
    end
    methods
        % Constructor
        function out = Tree(y,X,Leafmin,gamma,beta,temp)
            % Store data
            % out.X = X;
            % Create root node
            rootnode = Nodes(0,[],[],[],[],1:size(X,1),0);
            out.Allnodes{1} = rootnode;
            if ~isempty(Leafmin) && Leafmin > 0
                out.Leafmin = Leafmin;
            else
                out.Leafmin = 25;
            end
            out.NodeIds = 0;
            out.Smallnodes = 0;
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
            if gamma > 0 && gamma < 1
                out.gamma = gamma;
            else
                error('gamma must be between 0 and 1')
            end
            if beta >= 0
                out.beta = beta;
            else
                error('beta must be >= 0')
            end
            
            
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
            
            [~,out] = prior_eval(out,X);
            out.Temp = temp;
            out.Ntermnodes = 1;
            %out.Updateprior = 0;
        end
        
        % Calculate log-likelihood of a node and updates the node
        %   in the tree.
        function out = llike(obj,nodename,y)
            nind = nodeind(obj,nodename);
            out = obj;
            thenode = obj.Allnodes{nind};
            thenode = loglikefunc(thenode,obj,y);
            thenode.Updatellike = 0; 
            %out.Allnodes{nind}
            out.Allnodes{nind} = thenode;
        end
        
        % Calculate the likelihood of all the terminal
        %   nodes and update the tree log-likelihood
        function out = llike_termnodes(obj,y)
            out = obj;
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
        end
        
            
        % Birth Function
        function [out,birthindex] = birth(obj,y,X)
            % obj: an object of class "Tree"
            % y: a vector of the response variable
            % X: a matrix of the predictors
            
            out = obj;
            birthindex = [];
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

                node = obj.Allnodes{birthindex}; % birth node
                % Determine split rules, if necessary
                if node.Updatesplits == 1
                    node = getsplits(node,X,obj.Leafmin);
                end
                
                node = drawrule(node);
                
                if ~isempty(node.Rule)
                    % Put updated node in tree
                    out.Allnodes{birthindex} = node;  

                    % Add children nodes;
                    % Find unique new IDs for children
                    newids = newIDs(obj);
                    
                    if length(node.Xind) < out.Leafmin
                        error('bad start in leaf')
                    end
                    % Get data which is being passed down from parent
                    [XindL,XindR] = childrendata(out,birthID,X); % must use out since it has the new rule
                    
                    % TODO: remove error check
                    if length(XindL) < out.Leafmin || length(XindR) < out.Leafmin
                        error('leaf is bad')
                    end
                    
                    
                    % Add Nodes
                    newdepth = node.Depth + 1;
                    if any(ismember(newids,out.NodeIds))
                        error('New Ids are wrong.')
                    end
                    childnode1 = Nodes(newids(1),birthID,[],[],[],XindL,newdepth);
                    childnode2 = Nodes(newids(2),birthID,[],[],[],XindR,newdepth);
                    duplicateIDs(out);
                    out = addnode(out,childnode1,birthindex,'L',y);
                    out = addnode(out,childnode2,birthindex,'R',y);
                    duplicateIDs(out);
                    % Update log-likelihood
                    n = nnodes(out);
                    out.Lliketree = out.Lliketree - node.Llike + ...
                        out.Allnodes{n-1}.Llike + out.Allnodes{n}.Llike;
                    parentchildagree(out);
                    [~,out] = prior_eval(out,X);
                    out.Ntermnodes = out.Ntermnodes + 1;
                    return;
                end
            end
            warning('Birth Step Not Possible')
        end
        
        % Prune Function
        % pind is the pruned node's index in the original obj (tree)
        function [out,pind] = prune(obj,y,X)
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
                [~,out] = prior_eval(out,X);
                out.Ntermnodes = out.Ntermnodes - 1;
                parentchildagree(out);
                duplicateIDs(out);
            else
                error('Cannot prune a root node.')
            end
        end
        
        % Change Function (Currently with CART implementation - draw from prior)
        % p is the probabiliy of only changing the rule of a continuous
        % variable
        function [out,priordraw,startcont,endcont,nchange,nchange2] = change(obj,y,X,p)
            %cont = 0; % 1 if a continuous variable was changed
            nchange = []; % number of possible changes (up to 2) for continuous variables
            nchange2 = [];
            out0 = obj;
            % Find interior nodes
            [I,~] = interiornodes(obj);
            if isempty(I)
                error('Cannot perform a change step on a root node with no children.')
            end
            
            if length(I) > 1
                changeind_rand = randsample(I,length(I));
            else
                changeind_rand = I;
            end
            for ii = 1:length(I); % For each interior node
                changeind = changeind_rand(ii);
                node = out0.Allnodes{changeind};
                priordraw = 0;
                startcont = 0;
                endcont = 0;
                % Determine variable type of rule
                if strcmp(obj.Xclass{node.Rule{1}},'cell') % draw a rule from prior
                    priordraw = 1;
                else
                    startcont = 1;
                    if rand >= p
                        priordraw = 1; 
                    else
                        endcont = 1;
                    end
                end
                
                if 1
                    oldrule = node.Rule; % TODO: remove this later (error checking purpose)
                    % Update splits if necessary
                    if node.Updatesplits == 1
                        node = getsplits(node,X,obj.Leafmin);
                    end

                    % TODO: delete later...
                    % Error checking
                    psplits = node.Splitvals{oldrule{1}};
                    ruleinthere = 0;
                    if isa(psplits,'cell')
                        for jj = 1:length(psplits)
                            nrule = psplits{jj};
                            if length(nrule) == length(oldrule{2})
                                if all(strcmp(sort(nrule),sort(oldrule{2})))
                                    ruleinthere = 1;
                                end
                            end
                        end
                        if ~ruleinthere
                            
                            
                            
                            
                            myid = node.Id
                            rulevar = oldrule{1}
                            myrule = oldrule{2}
    %                         for jj = 1:length(out0.Allnodes)
    %                             jj
    %                             tnode = out0.Allnodes{jj};
    %                             [tnode.Id,tnode.Parent,tnode.Lchild,tnode.Rchild]
    %                             tnode.Rule
    %                             tnode.Rule{2}
    %                         end
                            printstructure(out0)
                            Treeplot(out0)
                            updatesplits = node.Updatesplits
                            % Print Rules of node
                            for mm = 1:length(node.Rule{2})
                                therules = node.Rule{2}{mm}
                            end

                            nodeupdated = getsplits(node,X,obj.Leafmin);

                            for mm = 1:length(nodeupdated.Rule{2})
                                therules2 = nodeupdated.Rule{2}{mm}
                            end

                            xindlength = length(node.Xind)
                            tabulate(X{node.Xind,node.Rule{1}})
                                    
                            % error('Old rule not a candidate for prior draw');
                        end
                    end
                end

                    
                
                
                % Store updated splitvals even if the change step does not
                % happen 
                out0.Allnodes{changeind} = node;
                
                if priordraw % Draw from the prior
                    % Choose a variable to split on
                    varind = find(node.nSplits > 0); % index on available variables
                    % Randomly order available variables
                    if length(varind) > 1
                        varind_rand = randsample(varind,length(varind));
                    elseif length(varind) == 1
                        varind_rand = varind;
                    else
                        error('No variables available (should not happen since original rule is valid)')
                    end
                    for jj = 1:length(varind_rand) % For each available variable
                        endcont = 0;
                        vind = varind_rand(jj);
                        rules = node.Splitvals{vind};
                        % Rearrange rules in random order
                        rules_rand = rules(randsample(length(rules),length(rules)));                            
                        for kk = 1:length(rules) % For each rule on this variable
                            if isa(rules_rand,'cell')
                                newrule = rules_rand{kk};
                            else
                                newrule = rules_rand(kk);
                                endcont = 1;
                            end
                            % For each possible rule
                            nodestar = node;
                            treestar = out0;
                            nodestar.Rule = {vind,newrule};
                            treestar.Allnodes{changeind} = nodestar;
                            % Check to see if rule leaves a tree with enough 
                            %   observations at each terminal node.
                            treestar = descendentdata(treestar,nodestar.Id,X);
                            
                            if treestar.Smallnodes == 0
                                % Compute log-likelihood of terminal nodes;
                                %treestar = llike(treestar,nodestar.Lchild,y);
                                %treestar = llike(treestar,nodestar.Rchild,y);
                                treestar = llike_termnodes(treestar,y);
                                out = treestar;
                                [~,out] = prior_eval(out,X);
                                % Fix any classification rules which are
                                % not up to date
                                for mm = 1:length(out.Allnodes)
                                    mynode = out.Allnodes{mm};
                                    if ~isempty(mynode.Rule)
                                        myrule = mynode.Rule{2};
                                        if isa(myrule,'cell')
                                            uq = unique(X{mynode.Xind,mynode.Rule{1}});
                                            myrule = myrule(ismember(myrule,uq)); % New Rule L
                                            mynode.Rule{2} = myrule;
                                            out.Allnodes{mm} = mynode;
                                        end
                                    end
                                end
                                parentchildagree(out);
                                duplicateIDs(out);
                                return;
                            end % otherwise try again
                        end
                    end
                else % Sequentially move rule one value to left or right
                    % Find the index of the rule in Splitvals
                    % Calcualte log-likelihood of chosen rule...
                    %cont = 1;
                    endcont = 1;
                    [n,trees] = getchangerules(out0,changeind,X);
                    nchange = sum(n);
                    if nchange == 2
                        if rand < .5
                            out = trees{1};
                        else
                            out = trees{2};
                        end
                    elseif nchange == 1
                        out = trees{n > 0};
                    else % Keep the same tree/rule;
                        out = out0;
                    end
                    out = llike_termnodes(out,y);
                    [~,out] = prior_eval(out,X);
                    n2 = getchangerules(out,changeind,X);
                    nchange2 = sum(n2);
                    return; % Never need to go through the loop
                end
                error('should never get here...');
            end
            error('No split value found. This should not happen since the original rule is valid.')
        end
        
        
        function [n,trees] = getchangerules(obj,nodeind,X)
            trees = cell(2,1);
            n = [0,0];
            node = obj.Allnodes{nodeind};
            vind = node.Rule{1};
            rules = node.Splitvals{vind};
            nrules = length(rules);
            oldrule = node.Rule{2};
            rind = find(oldrule == rules);
            if isempty(rind) % possibly rules have changed above it...
                % This does not change the leaf nodes below...
                rind = find(oldrule >= rules,1,'last');
            end
            Ltreecalc = 0;
            Rtreecalc = 0;
            if rind > 1 && rind < nrules; % If rule isn't on boundary
                Ltreecalc = 1;
                Rtreecalc = 2;
            elseif rind == 1 && nrules > 1
                Rtreecalc = 1;
            elseif rind == nrules && nrules > 1
                Ltreecalc = 1;
            elseif rind == 1 && nrules == 1
                % Move either direction isn't possible
            else
                rind
                nrules
                error('Logic is bad.');
            end
            if Ltreecalc
                newruleL = {vind,rules(rind - 1)};
                treeL = obj;
                nodeL = node;
                nodeL.Rule = newruleL;
                treeL.Allnodes{nodeind} = nodeL;
                treeL = descendentdata(treeL,nodeL.Id,X);
                if treeL.Smallnodes == 0
                    n(1) =  1;
                    trees{1} = treeL;
                end 
            end
            if Rtreecalc
                newruleR = {vind,rules(rind + 1)};
                treeR = obj;
                nodeR = node;
                nodeR.Rule = newruleR;
                treeR.Allnodes{nodeind} = nodeR;
                treeR = descendentdata(treeR,nodeR.Id,X);
                if treeR.Smallnodes == 0
                    n(2) = 1;
                    trees{2} = treeR;
                end 
            end
        end
                        
        
        
        % Swap
        % If childID is specified, a swap is done with the childID and it's parent. 
        function [out, swappossible] = swap(obj,y,X,childID,llikepriorcalc)
            if isempty(childID)
                % Find internal nodes whose parent is also a internal node;
                Ids = parentchildpairs(obj);
                if isempty(Ids)
                    error('Cannot perform swap step: No internal parent-child pairs.')
                end
                % Randomize Ids
                Ids = Ids(randsample(length(Ids),length(Ids)));
                swappossible = [];
            else
                if length(childID) > 1
                    error('childID must be of length 1');
                end
                childIDind = nodeind(obj,childID);
                tnode = obj.Allnodes{childIDind};
                if isempty(tnode.Parent) || isempty(tnode.Rule)
                    error('Child node must be an interior node with a parent');
                end
                Ids = childID;
                swappossible = 0;
            end
            
            % Loop through until we find a swap without empty (or too small of) nodes.
            for ii = 1:length(Ids)
                out = obj;
                % Randomly select one node
                rID = Ids(ii);
                nind = nodeind(out,rID);
                node = out.Allnodes{nind};
                % Child's rule
                crule = node.Rule; % child rule
                % Find the node's parent
                nindParent = nodeind(out,node.Parent);
                nodeParent = out.Allnodes{nindParent};
                % Parent Rule
                prule = nodeParent.Rule;
                % Find grandparent
                grandparent = 0;
                if ~isempty(nodeParent.Parent)
                    nindGrandparent = nodeind(out,nodeParent.Parent);
                    nodeGrandparent = out.Allnodes{nindGrandparent};
                    grandparent = 1;
                end
                % Do the parent and child have a rule on the same variable?
                if prule{1} ~= crule{1} || ((prule{1} == crule{1}) &&  strcmp(out.Xclass(prule{1}),'cell'))% If no, Swap.
                    % Check if both children have same split rule
                    % Is the node a left or right child?
                    LRind = ismember([nodeParent.Lchild,nodeParent.Rchild],rID);
                    indR = nodeind(out,nodeParent.Rchild);
                    indL = nodeind(out,nodeParent.Lchild);
                    if all(LRind == [1 0]) % Left Child
                        lrule = node.Rule;
                        %indR = nodeind(obj,nodeParent.Rchild);
                        rrule = obj.Allnodes{indR}.Rule;
                        ind2 = indR;
                    elseif all(LRind == [0 1]) % Right Child
                        rrule = node.Rule;
                        %indL = nodeind(obj,nodeParent.Lchild);
                        lrule = obj.Allnodes{indL}.Rule;
                        ind2 = indL;
                    end

                    % Are the two children rules equal?
                    eqrule = 0;
                    if ~isempty(lrule) && ~isempty(rrule)
                        if lrule{1} == rrule{1}
                            if isa(lrule{2},'cell')
                                if isequal(sort(lrule{2}),sort(rrule{2}))
                                    eqrule = 1;
                                end
                            else % Numeric
                                if lrule{2} == rrule{2}
                                    eqrule = 1;
                                end
                            end                            
                        end
                    end
                    % If not equal rules, just swap one
                    % Swap rules
                    nodeParent.Rule = crule;
                    node.Rule = prule;
                    % Put rules in tree
                    out.Allnodes{nind} = node;
                    out.Allnodes{nindParent} = nodeParent;
                    % change rule of other child if equal rules (CART)
                    if eqrule
                        node2 = out.Allnodes{ind2};
                        % Update rule
                        node2.Rule = prule;
                        % Update tree
                        out.Allnodes{ind2} = node2;
                    end    
                     % Update descendent data
                    out = descendentdata(out,nodeParent.Id,X);
                    if out.Smallnodes == 0 % we have enough data at each terminal node
                        if llikepriorcalc
                            out = llike_termnodes(out,y);
                            [~,out] = prior_eval(out,X);
                        end
                        % Fix any classification rules which are
                        % not up to date
                        for mm = 1:length(out.Allnodes)
                            mynode = out.Allnodes{mm};
                            if ~isempty(mynode.Rule)
                                myrule = mynode.Rule{2};
                                if isa(myrule,'cell')
                                    uq = unique(X{mynode.Xind,mynode.Rule{1}});
                                    myrule = myrule(ismember(myrule,uq)); % New Rule L
                                    mynode.Rule{2} = myrule;
                                    out.Allnodes{mm} = mynode;
                                end
                            end
                        end
                        parentchildagree(out);
                        duplicateIDs(out);
                        swappossible = 1;
                        return;
                    end 
                else % Rotate if rule on same variable.
                    % TODO: Add in the rotate steps. (most complex).
                    % Is the node a left or right child?
                    LRind = ismember([nodeParent.Lchild,nodeParent.Rchild],rID);
                    indR = nodeind(out,nodeParent.Rchild);
                    indL = nodeind(out,nodeParent.Lchild);
                    if all(LRind == [1 0]) % Left Child
                        %lrule = node.Rule;
                        %indR = nodeind(obj,nodeParent.Rchild);
                        %rrule = obj.Allnodes{indR}.Rule;
                        LR = 'L';
                    elseif all(LRind == [0 1]) % Right Child
                        %rrule = node.Rule;
                        %indL = nodeind(obj,nodeParent.Lchild);
                        %lrule = obj.Allnodes{indL}.Rule;
                        LR = 'R';
                    end
                    if strcmp(LR,'L') % Do a right rotation
                        % Original parents/children
                        pparent = nodeParent.Parent;
                        crchild = node.Rchild;
                        nind2 = nodeind(out,crchild);
                        node2 = out.Allnodes{nind2};
                        % Change parents and children
                        nodeParent.Parent = node.Id; 
                        nodeParent.Lchild = crchild;
                        node.Parent = pparent;
                        node.Rchild = nodeParent.Id;
                        node2.Parent = nodeParent.Id;
                    else % Do a left rotation
                        pparent = nodeParent.Parent;
                        clchild = node.Lchild;
                        nind2 = nodeind(out,clchild);
                        node2 = out.Allnodes{nind2};
                        % Change parents and children
                        nodeParent.Parent = node.Id;
                        nodeParent.Rchild = clchild;
                        node.Parent = pparent;
                        node.Lchild = nodeParent.Id;
                        node2.Parent = nodeParent.Id;        
                    end
                    if grandparent
                        if nodeGrandparent.Lchild == nodeParent.Id
                            nodeGrandparent.Lchild = node.Id;
                        elseif nodeGrandparent.Rchild == nodeParent.Id
                            nodeGrandparent.Rchild = node.Id;
                        else
                            error('Logic is bad.')
                        end
                        % Put grandparent node back in tree
                        out.Allnodes{nindGrandparent} = nodeGrandparent;
                    end
                    nodeParent.Depth = nodeParent.Depth + 1;
                    node.Depth = node.Depth - 1;
                    % Assign the Xind from parent to child
                    node.Xind = nodeParent.Xind;
                    
                    % The splits will need to be updated
                    node.Updatesplits = 1;
                    nodeParent.Updatesplits = 1;
                    node.Updatellike = 1;
                    nodeParent.Updatellike = 1;

                    % Put nodes back in tree
                    out.Allnodes{nind} = node;
                    out.Allnodes{nindParent} = nodeParent;
                    out.Allnodes{nind2} = node2;
                    % Update descendent data
                    if grandparent
                        out = descendentdata(out,nodeGrandparent.Id,X);
                    else
                        node.Xind = 1:size(X,2);
                        out = descendentdata(out,node.Id,X);
                    end
                    % Update Depths
                    out = depthupdate(out,node.Id);
                    parentchildagree(out);
                    duplicateIDs(out);
                    swappossible = 1;
                    if llikepriorcalc
                        [~,out] = prior_eval(out,X);
                    end
                    return;
                end
            end  
            out = obj; % Return same tree
            % warning('Swap step not possible.')
        end
        
        % Update the depths of all nodes below node with Id
        function out = depthupdate(obj,Id)
            out = obj;
            nind = nodeind(out,Id);
            node = out.Allnodes{nind};
            id_L = node.Lchild;
            id_R = node.Rchild;
            nindL = nodeind(out,id_L);
            nindR = nodeind(out,id_R);
            nodeL = out.Allnodes{nindL};
            nodeR = out.Allnodes{nindR};
            nodeL.Depth = node.Depth + 1;
            nodeR.Depth = node.Depth + 1;
            out.Allnodes{nindL} = nodeL;
            out.Allnodes{nindR} = nodeR;
            if ~isempty(nodeL.Rule)
                out = depthupdate(out,id_L);
            end
            if ~isempty(nodeR.Rule)
                out = depthupdate(out,id_R);
            end
        end
            
        
        
%         % Swap
%         function out = swap(obj,y,X)
%             % Find internal nodes whose parent is also an internal node;
%             Ids = parentchildpairs(obj);
%             if isempty(Ids)
%                 error('Cannot perform swap step: No internal parent-child pairs.')
%             end
%             Ids = randsample(Ids,length(Ids));
%             % Loop through until we find a swap without empty nodes.
%             for ii = 1:length(Ids)
%                 tmpout = obj;
%                 % Randomly select one child
%                 rind = Ids(ii);
%                 nind = nodeind(obj,rind);
%                 node = obj.Allnodes{nind};
%                 % Find the node's parent
%                 nindParent = nodeind(obj,node.Parent);
%                 nodeParent = obj.Allnodes{nindParent};
%                 % Parent Rule
%                 prule = nodeParent.Rule;
%                 % Is the node a left or right child?
%                 LRind = ismember([nodeParent.Lchild,nodeParent.Rchild],rind);
%                 indR = nodeind(obj,nodeParent.Rchild);
%                 indL = nodeind(obj,nodeParent.Lchild);
%                 if all(LRind == [1 0]) % Left Child
%                     lrule = node.Rule;
%                     %indR = nodeind(obj,nodeParent.Rchild);
%                     rrule = obj.Allnodes{indR}.Rule;
%                 elseif all(LRind == [0 1]) % Right Child
%                     rrule = node.Rule;
%                     %indL = nodeind(obj,nodeParent.Lchild);
%                     lrule = obj.Allnodes{indL}.Rule;
%                 end
% 
%                 % Are the two children rules equal?
%                 eqrule = SplitRule.equalrule(lrule,rrule);
%                 if eqrule
%                     % Switch Rules
%                     newprule = lrule;
%                     newLrule = prule;
%                     newRrule = prule;
%                 else
%                     if all(LRind == [1 0]) % Left Child
%                         newprule = lrule;
%                         newLrule = prule;
%                         newRrule = rrule;
%                     elseif all(LRind == [0 1]) % Right Child
%                         newprule = rrule;
%                         newLrule = lrule;
%                         newRrule = prule;
%                     end
%                 end
%                 tmpout.Allnodes{nindParent}.Rule = newprule;
%                 tmpout.Allnodes{indL}.Rule = newLrule;
%                 tmpout.Allnodes{indR}.Rule = newRrule;
% 
%                 % Update Descendent Data
%                 tmpout.Smallnodes = 0;
%                 tmpout = descendentdata(tmpout,nodeParent.Id,X);
%                 if tmpout.Smallnodes == 0
%                     out = tmpout;
%                     % Update log-likelihood on terminal nodes which need updating
%                     % Also update Lliketree
%                     out.Lliketree = 0;
%                     [I,Ids] = termnodes(out);
%                     for jj=1:length(Ids)
%                         if out.Allnodes{I(jj)}.Updatellike
%                              out = llike(out,Ids(jj),y);
%                         end
%                         out.Lliketree = out.Lliketree + out.Allnodes{I(jj)}.Llike;
%                     end
%                     return;
%                 end
%             end
%             error('Swap step not possible.')
%             
%         end
            
              
        % Find interior nodes (Ids) with a parent who is also an interior node;
        function Ids = parentchildpairs(obj)
            [~,Id] = termnodes(obj);
            nIds = obj.NodeIds;
            %nIds = nIds(nIds > 0); % Exclude the root node;
            for ii=1:length(nIds)
                if isempty(obj.Allnodes{ii}.Parent)
                    rootid = obj.Allnodes{ii}.Id;
                end
            end
            nIds = nIds(nIds ~= rootid);
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
           node.Rule = [];
           
           newind = nodeind(out,node.Id);
           out.Allnodes{newind} = node; % update node
           % Update likelihood if necessary
           if node.Updatellike
               out = llike(out,out.Allnodes{newind}.Id,y);
           end
           % Add in Log-likelihood of the parent node
           out.Lliketree = out.Lliketree + out.Allnodes{newind}.Llike;
        end
        
        % Add node to tree, and calculate log-likelhood of node
        function out = addnode(obj,node,parentind,LR,y)
            n = nnodes(obj);
            out = obj;
            out.Allnodes{n + 1} = node;
            out.NodeIds(n + 1) = node.Id;
            if strcmp(LR,'L')
                out.Allnodes{parentind}.Lchild = node.Id;
            elseif strcmp(LR,'R')
                out.Allnodes{parentind}.Rchild = node.Id;
            else
                error('LR arg must be either "L" or "R"');
            end
            % Calculate log-likelihood of new node
            out = llike(out,node.Id,y);
        end
        
        % Returns the index and Ids of the interior nodes
        function [I,Id] = interiornodes(obj)
            n = nnodes(obj);
            if n <= 1
               I = [];
               Id = [];
               return;
            end
            Id = zeros(n,1);
            I = Id;
            cntr = 1;
            for ii=1:n
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
            n = nnodes(obj);
            Id = zeros(n,1);
            I = Id;
            cntr = 1;
            for ii=1:n
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
            %n = nnodes(obj);
            %allIds = zeros(n,1);
            %for ii=1:n
            %    allIds(n) = obj.Allnodes{ii}.Id;
            %end
            allIds = obj.NodeIds;
            %alln = 0:(n-1);
            alln = 0:max(allIds);
            I = ~ismember(alln,allIds);
            out = alln(I);
            cntr = 0;
            n = max(allIds) + 1;
            while length(out) < 2
                out = [out,n + cntr];
                cntr = cntr + 1;
            end
        end
        
        % Returns the number of nodes (including leaves and root)
        function out = nnodes(obj)
           out = length(obj.Allnodes); 
        end
        
        % Obtain the data for each of the children of a designated node
        % ndata is a vector of length 2 giving the number of data points in
        % the left and right child nodes, respectively.
        function [XindL,XindR,ndata] = childrendata(obj,nodeid,X)
            nind = nodeind(obj,nodeid);
            node = obj.Allnodes{nind};
            therule = node.Rule; 
            if isempty(therule)
                error('No rule specified for the node.')
            end
            if strcmp(obj.Xclass(therule{1}),'double')
                %I = find(table2array(X(node.Xind,therule.Varcol)) <= therule.Varrule);
                % I = table2array(X(node.Xind,therule{1})) <= therule{2} ;
                I = X{node.Xind,therule{1}} <= therule{2};
            elseif strcmp(obj.Xclass(therule{1}),'cell')
                %I = ismember(table2cell(X(node.Xind,therule{1})),therule{2}{1});
                % I = ismember(table2cell(X(node.Xind,therule{1})),therule{2});
                I = ismember(X{node.Xind,therule{1}},therule{2});
            else
                error('Unexpected variable class found.')
            end
            XindL = node.Xind(I);
            XindR = node.Xind(~I); 
            ndata = [sum(I),sum(~I)];
            
            % Check for children with data counts less than leafmin
            
%             
%             % Check for empty children (children with no data)
%             emptyL = 0;
%             emptyR = 0;
%             if min(size(XindL)) == 0
%                 emptyL = 1;
%                 XindL = [];
%             end
%             if min(size(XindR)) == 0
%                 emptyR = 1;
%                 XindR = [];
%             end
%             empty = [emptyL,emptyR];               
            
        end
        
        % Update data (Xind) of all descendents of a node (recursive function)
        function out = descendentdata(obj,nodeid,X)
            nind = nodeind(obj,nodeid);
            node2 = obj.Allnodes{nind};
            if ~isempty(node2.Lchild) && ~isempty(node2.Rchild)
                [XindL,XindR,ndata] = childrendata(obj,nodeid,X);           
                out = obj;
                if any(ndata < obj.Leafmin)
                    out.Smallnodes = 1;
                end
                nindL = nodeind(out,node2.Lchild);
                nindR = nodeind(out,node2.Rchild);
                out = updatedata(out,nindL,XindL);
                out = updatedata(out,nindR,XindR);
                
                
                
%                 nodeL = out.Allnodes{nindL};
%                 nodeR = out.Allnodes{nindR};
%                 if ~isempty(nodeL.Rule)
%                     if isa(nodeL.Rule{2},'cell')
%                         Treeplot(out)
%                         ruleL = nodeL.Rule{2}
%                         uq = unique(X{nodeL.Xind,nodeL.Rule{1}})
%                         ruleL = ruleL{ismember(ruleL,uq)}; % New Rule L
%                         nodeL.Rule{2} = ruleL;
%                         out.Allnodes{nindL} = nodeL;
%                     end
%                 end
%                 if ~isempty(nodeR.Rule)
%                     if isa(nodeR.Rule{2},'cell')
%                         ruleL = nodeR.Rule{2};
%                         uq = unique(X{nodeR.Xind,nodeR.Rule{1}});
%                         ruleR = ruleR{ismember(ruleR,uq)}; % New Rule R
%                         nodeR.Rule{2} = ruleR;
%                         out.Allnodes{nindR} = nodeR;
%                     end
%                 end
                
                
                
                % Now do it for the descendents
                out = descendentdata(out,node2.Lchild,X);
                out = descendentdata(out,node2.Rchild,X);
            else
                out = obj;
            end
        end
        
        
        
        % Get the index of a node with nodeid;
        function out = nodeind(obj,nodeid)
            n = length(nodeid);
            if n == 1
                out = find(obj.NodeIds == nodeid);
            elseif n > 1
                out = zeros(n,1);
                for ii=1:n
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
            % If only the root node, return empty things
            if length(I) < 2 && length(obj.Allnodes) < 2
                %I = 1;
                %ids = obj.Allnodes{I}.Id;
                I = [];
                ids = [];
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
                node1.Updatesplits = 1;
            elseif ~all(sort(Xind) == sort(node1.Xind))
                node1.Updatellike = 1;
                node1.Updatesplits = 1;
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
        
        % Evaluate the prior on the tree (obj)
        %  Will also return the tree with updated split values
        function [lprior,out] = prior_eval(obj,X)
            out = obj;
            lprior = 0;
            for ii=1:length(obj.Allnodes)
                % prior on depth
                node = obj.Allnodes{ii};
                d = node.Depth;
                if ~isempty(node.Rule) % if not a terminal node
                    if node.Updatesplits % update splits if needed
                        node = getsplits(node,X,obj.Leafmin);
                        out.Allnodes{ii} = node;
                    end
                    lprior = lprior + log(obj.gamma) - obj.beta*log(1 + d) - ...
                        log(sum(node.nSplits > 0)) - log(node.nSplits(node.Rule{1}));
                    if ~isfinite(lprior)
                        betadepth = obj.beta*log(1 + d)
                        nvars = log(sum(node.nSplits > 0))
                        nsplits = log(node.nSplits(node.Rule{1}))
                        node.nSplits(node.Rule{1})
                        node.Splitvals{node.Rule{1}}
                        Treeplot(out)
                        id = node.Id
                        therule = node.Rule{2}
                        printstructure(out)
                        error('non-finite prior evaluation')
                    end
                %log(sum(node.nSplits));
                else % if a terminal node
                    lprior = lprior + log(1 - obj.gamma/(1 + d)^obj.beta);
                end                   
            end
            out.Prior = lprior;
        end
            
        
        
        % Graphs
        function treelines(obj,nodename,level,treedepth,parentxloc,LR,plotdens,y)
            width = 1; % space between terminal nodes
            % maxloc = 1 + width*(2*treedepth - 1);
            nind = nodeind(obj,nodename);
            node = obj.Allnodes{nind};
            if ~isempty(plotdens)
                adddens = 1;
            else
                adddens = 0;
            end
            if ~isempty(node.Parent)
            %if nodename > 0
                if LR == 'L'
                    plusminus = -1;
                elseif LR == 'R'
                    plusminus = 1;
                else
                    error('Must correctly specify "LR" as "L" or "R"')
                end
                %delta = width*2^-(treedepth + level - 1);
                delta = width*2^(treedepth + level - 1); % Old delta
                % which works...
                %delta = width - 1/(2*treedepth);
                %delta = (1 + width*(2*treedepth - 1) - (maxloc+1)/2)/treedepth;
                xval = parentxloc + plusminus*delta;
                if ~adddens
                    plot([parentxloc, xval],[level + 1,level],'k')
                end
            else
                % Get the xval to pass to children
                %maxloc = 2*treedepth*width;     
                % xval = (maxloc + 1)/2;
                %parentxloc = xval; % just for the root node...
                xval = 0;
            end
            
            % Plot the rule of the parent
            if ~isempty(node.Rule) % if parent node
                colnum = node.Rule{1};
                colname = obj.Varnames(colnum);
                if strcmp(obj.Xclass(node.Rule{1}),'double')
                    ruletext = strcat(colname,' \leq ',num2str(node.Rule{2}));
                else
                    for ii=1:length(node.Rule{2})
                        if ii == 1
                            grp = node.Rule{2}{ii};
                        else
                            %grp = strcat(grp,' , ',node.Rule{2}{ii});
                            grp = [grp,', ',node.Rule{2}{ii}];
                        end
                    end
                    ruletext = strcat(colname,' \in \{',grp,'\}');
                    %ruletext = strcat(colname,' \in \{',cell2mat(node.Rule.Varrule),'\}');
                end
                if ~adddens
                    text(xval,level,ruletext,'HorizontalAlignment','right')
                end
            end
            
            %n = length(obj.Allnodes);
            if ~isempty(node.Lchild) && ~isempty(node.Rchild)
                treelines(obj,node.Lchild,level-1,treedepth,xval,'L',plotdens,y)
                treelines(obj,node.Rchild,level-1,treedepth,xval,'R',plotdens,y)
            elseif isempty(node.Lchild) && isempty(node.Rchild) && ~isempty(plotdens)
                %hold off
                % Now plot the density at the terinal nodes
                % delta2 = width*2^(treedepth + level);
                %delta = width - 1/(2*treedepth);
                %delta = (1 + width*(2*treedepth - 1) - (maxloc+1)/2)/treedepth;
                %axisparms = gca;
                %xrange = axisparms.XLim(2) - axisparms.XLim(1);
                %yrange = axisparms.YLim(2) - axisparms.YLim(1);
                xrange = plotdens(1,2) - plotdens(1,1);
                yrange = plotdens(2,2) - plotdens(2,1);
                xorigin = plotdens(1,1);
                yorigin = plotdens(2,1);
                xvald = (xval - xorigin)/xrange;
                yvald = (level - yorigin)/yrange;
                % Normalize again
                lw = 2^(-treedepth)*plotdens(4,1);
                xvald = xvald*plotdens(4,1) + plotdens(3,1) - lw/2;
                yvald = yvald*plotdens(4,2) + plotdens(3,2) - lw;
                % [xvald,yvald]
                axes('Position',[xvald,yvald,lw,lw])
                box on
                if isempty(y)
                    error('No data to plot densities.')
                else
                    %size(node.Xind)
                    %[node.Id, nind]
                    try 
                        lgpdens(y(node.Xind))
                    catch
                        lgpdens(y(node.Xind),'imp_sampling','off')
                    end
                    
                end
                %hold on
            else
                % error('Unexpected case encountered.')
            end
        end
        
        function Treeplot(obj,y)
            n = length(obj.Allnodes);
            if n <= 1
                error('Tree must be more than a root node.')
            else
                nodegraph = zeros(1,n);
                for ii = 1:n
                    parent = obj.Allnodes{ii}.Parent;
                    if ~isempty(parent)
                        % nodegraph(ii) = parent + 1;
                        nodegraph(ii) = obj.Allnodes{ii}.Depth;
                    else
                        rootnodename = obj.Allnodes{ii}.Id;
                    end
                end
                %treeplot(nodegraph)
                treedepth = max(nodegraph);
                %width = 1; % space between terminal nodes
                %maxloc = 2*treedepth;
                figure;
                hold on;
                treelines(obj,rootnodename,0,treedepth,'','',[],[])
                %hold off;
                axisparms = gca;
                plotdens = [axisparms.XLim; axisparms.YLim;
                    axisparms.Position(1:2); axisparms.Position(3:4)];
                axis off;
                treelines(obj,rootnodename,0,treedepth,'','',plotdens,y)
                hold off;
            end
        end
        
        % Troubleshooting functions
        % Prind ID and ID status
        function llikestatus(obj)
            for ii = 1:length(obj.Allnodes)
                node = obj.Allnodes{ii};
                if isempty(node.Rule)
                    tnode = 'yes';
                else
                    tnode = 'no';
                end
                disp(['Id=',num2str(node.Id),', Updatellike=',...
                    num2str(node.Updatellike),...
                    ', Updatesplits=',num2str(node.Updatesplits),...
                    ', Terminalnode=',tnode]);
            end     
        end
        
        function printstructure(obj)
            for ii = 1:length(obj.Allnodes)
                node = obj.Allnodes{ii};
                if isempty(node.Rule)
                    tnode = 'yes';
                else
                    tnode = 'no';
                end
                if ~isempty(node.Rule)
                    therule = node.Rule{2};
                    if isa(therule,'double')
                        therule = num2str(therule);
                    else
                        therule = char(therule)';
                    end
                    rulevar = num2str(node.Rule{1});
                else
                    rulevar = '';
                    therule = '';
                end
                disp(['Id=',num2str(node.Id),...
                    ', Parent=',num2str(node.Parent),...
                    ', Lchild=',num2str(node.Lchild),...
                    ', Rchild=',num2str(node.Rchild),...
                    ', Rulevar=',rulevar,...
                    ', Rule=',therule,...
                    ', Terminalnode=',tnode]);
                fprintf('\n')
            end     
        end
        
        function parentchildagree(obj)
            for ii = 1:length(obj.Allnodes)
                node = obj.Allnodes{ii};
                if ~isempty(node.Lchild)
                    nindL = nodeind(obj,node.Lchild);
                    if obj.Allnodes{nindL}.Parent ~= node.Id
                        etext = ['Node ',num2str(node.Id),...
                            ' left child does not match.'];
                        error(etext)
                    end
                end
                if ~isempty(node.Rchild)
                    nindR = nodeind(obj,node.Rchild);
                    if obj.Allnodes{nindR}.Parent ~= node.Id
                        etext = ['Node ',num2str(node.Id),...
                            ' right child does not match.'];
                        error(etext)
                    end
                end
            end
        end
        
        function duplicateIDs(obj)
            if length(obj.Allnodes) ~= length(unique(obj.NodeIds))
                error('Duplicate Ids encountered.')
            end
        end
        
        function [n,out] = nbirthnodes(obj,X)
            out = obj;
            n = 0;
            for ii = 1:length(obj.Allnodes)
                if isempty(obj.Allnodes{ii}.Rule) % terminal node
                    node = obj.Allnodes{ii};
                    if node.Updatesplits
                        node = getsplits(node,X,obj.Leafmin);
                        out.Allnodes{ii} = node;
                    end
                    if sum(node.nSplits) > 0
                        n = n + 1;
                    end
                end
            end
        end
        
        function [nnodes,out] = nchangenodes(obj,X)
            out = obj;
            nnodes = 0; % Total number of change nodes
            for ii = 1:length(obj.Allnodes) % For all nodes
                if ~isempty(obj.Allnodes{ii}.Rule) % If not a terminal node
                    node = obj.Allnodes{ii};
                    if node.Updatesplits
                        node = getsplits(node,X,obj.Leafmin);
                        out.Allnodes{ii} = node;
                    end
                    if sum(node.nSplits) > 0
                        nnodes = nnodes + 1;
                    end
                end
            end         
        end
        % Count the number of nodes which can be swapped
        function n = nswaps(obj,y,X)
            n=0;
            for ii = 1:length(obj.Allnodes)
                node = obj.Allnodes{ii};
                if ~isempty(node.Rule) && ~isempty(node.Parent) % swapable
                    [~,sp] = swap(obj,y,X,node.Id,0);
                    n = n + sp;
                end
            end
        end
        
        % Check to see if the getsplits funcion has failed
        function [problem,I] = getsplitcheck(obj,X)
            problem = 0;
            I = [];
            for ii = 1:length(obj.Allnodes)
                node = obj.Allnodes{ii};
                updatednode = getsplits(node,X,obj.Leafmin);
                if node.Updatesplits == 0
                    if updatednode.Updatesplits == 0
                        % Check to make sure they match
                        if all(node.nSplits == updatednode.nSplits)
                            for jj = 1:length(node.Splitvals)
                                rules1 = sort(node.Splitvals{jj});
                                rules2 = sort(updatednode.Splitvals{jj});
                                if isa(rules1,'cell')
                                    if ~all(strcmp(rules1,rules2))
                                        problem = 1;
                                        I = [I,ii];
                                    end
                                else
                                    if ~all(rules1 == rules2)
                                        problem = 1;
                                        I = [I,ii];
                                    end
                                end   
                            end
                        else
                            problem = 1;
                            I = [I,ii];
                        end
                    else
                        error('Should not happen.')
                    end
                end
            end
        end
            
        
    end             
end