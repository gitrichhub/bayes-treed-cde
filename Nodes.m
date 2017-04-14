classdef Nodes
    properties
        Depth % How deep is the node in the tree?
        Id % Id for the node
        Parent % Id for the parent node
        Lchild % Id of left child
        Rchild % Id of right child
        Rule % Split rule for the node, cell of length 2 which
             %    contains the variable number/column and rule.
        % AvailableSplits % vector of size p (# of predictors) by 1 
                       %    indicating the number of possible
                       %    splits.
        Xind % Index of which data values pass through it.
        Llike % Log-Likelihood contribution of data in the partition...
        Splitvals
        nSplits
        Updatellike % 0 if no llike update/computation needed, 
                       % 1 if it is needed
        Updatesplits % 0 if no update of the splits is necessary, 1 if needed
    end
    methods
        % Constructor
        function obj = Nodes(Id,Parent,Lchild,Rchild,Rule,Xind,Depth) % Add in an Xind later?
            % Id: ID of the new node (integer)
            % Parent: ID of parent (integer)
            % Lchild: ID of left child (integer)
            % Rchild: ID of right child (integer)
            % Rule: A SplitRule class object defining the rule for the node
            % Xind: An index of which data points 
            % Depth: How many levels deep the rule is 
            %          (root node is depth 1)
            
            % Assign ID
            if isa(Id,'numeric')
                obj.Id = Id;
            else
                error('Id must be numeric')
            end
            
            % Assign Parent
            if isa(Parent,'numeric')
                obj.Parent = Parent;
            else
                error('Parent Id must be numeric')
            end
            
            % Assign Children
            if isa(Lchild,'numeric')
                obj.Lchild = Lchild;
            else
                error('Lchild must be numeric')
            end
            if isa(Rchild,'numeric')
                obj.Rchild = Rchild;
            else
                error('Rchild must be numeric')
            end
            
            % Assign the split rule for the node
            if iscell(Rule)
                if length(Rule) == 2
                    obj.Rule = Rule;
                else 
                    error('Rule must be of length 2.')
                end
            elseif isempty(Rule)
            else 
                error('Rule must be empty or a cell object of length 3.')
            end
            
            % Assign the index for which data values pass through this
            %   node.
            if isa(Xind,'numeric')
                obj.Xind = Xind;
            elseif isempty(Xind)
            else
                error('Xind must be numeric or empty')
            end
            
%             if isa(Leafmin,'numeric')
%                 if isempty(Leafmin)
%                     Leafmin = 25;
%                 elseif Leafmin > 0
%                 else
%                     error('Leafmin must be greater than 0')
%                 end
%             else
%                 error('Leafmin must be an integer greater than 0')
%             end
            
            % Assign Available Predictor Matrix
%             if isa(X,'table')
%                 % Nodes.createpredmatrix(X,Xind,Leafmin);
%                 % getSplits(X,Xind,Leafmin)
%             else
%                 error('X must be a table')
%             end
            
            % Assign Depth
            if isa(Depth,'numeric')
                obj.Depth = Depth;
            else
                error('Depth must be numeric')
            end
            
            % Llike needs to be computed since it hasn't been computed yet
            obj.Updatellike = 1;
            obj.Updatesplits = 1;
        end
        
        
        
        
        % Add/Change Rule
%         function out = newrule(obj,therule)
%             if isa(therule,'SplitRule') || isempty(therule)
%                 obj.Rule = therule;
%                 out = obj;
%             else
%                 error('rule must be of class "SplitRule"');
%             end
%         end
        
%         % select rule for factors
%         function out = factrule(obj)
%             out.Xind
%         end
        
        % Draw a new rule from the prior
        function out = drawrule(node) % obj is a Tree class object
            out = node;
            if out.Updatesplits == 1 % Update splits if necessary
                error('Splits must be updated prior to calling drawrule.')
            end
            
            if sum(out.nSplits) > 0
                vindex = find(out.nSplits > 0);
                vind = vindex(randsample(length(vindex),1)); % Variable index
                % Now randomly sample a rule
%                 nsplits = out.nSplits
%                 splitvals = out.Splitvals
%                 rule = out.Rule
%                 parent = out.Parent
%                 lchild = out.Lchild
%                 rchild = out.Rchild
                svals = out.Splitvals{vind};
                if isa(svals,'cell')
                    newrule = svals{randsample(length(svals),1)};
                else
                    newrule = svals(randsample(length(svals),1));
                end
                out.Rule = {vind,newrule};
            end % else leave rule empty
        end
        
        
        
        % Draw a new rule from the prior
%         function [nrule, colind] = drawrule(node,obj,X,colind) % obj is a Tree class object
%             % Choosing colind = 0 will randomly choose a variable to split
%             % Randomly choose a column to split on
%             if colind < 1
%                 % DO SOME CHECKING BEFORE CHOOSING THE COLUMN:
%                 %   check to see if there are at least two categories
%                 %   check to make sure we have enough data to justify the
%                 %   split.
%                 % Randomly choose a column to split on (if not specified)
%                 candcols = [];
%                 for ii = 1:size(X,2)
%                     %if strcmp(obj.Xclass(ii),'double')
%                         nu = max(size(unique(X(node.Xind,ii))));
%                         if nu > 1
%                             candcols = [candcols, ii];
%                         end
%                     %elseif strcmp(obj.Xclass(ii),'cell')
%                 end
%                 if ~isempty(candcols)
%                     colind = randsample(candcols,1);
%                 else
%                     nrule = [];
%                     colind = 0;
%                     return
%                 end
%             end
%             if strcmp(obj.Xclass(colind),'double')
%                 rulevals = randsample(table2array(X(node.Xind,colind)),length(node.Xind));
%                 
%                 % OLD VERSION...
%                 % Don't allow the rule to be the maximum (gives an empty
%                 % node)
%                 mval = max(rulevals);
%                 for ii = 1:length(rulevals)
%                     if rulevals(ii) < mval
%                         nrule = SplitRule(X,colind,rulevals(ii));
%                         break
%                     end
%                 end
%                 
% 
%                 % New Version
%                 % Subset rules so that you have at least leafmin obs at
%                 %   each terminal node...
% %                 rulevals = sort(table2array(X(node.Xind,colind)));
% %                 len = length(rulevals);
% %                 leafmin = obj.Leafmin;
% %                 if len >= 2*leafmin
% %                     rulevals = rulevals([leafmin,len-leafmin]);
% %                     if length(rulevals) > 1
% %                         ruleval = randsample(ruleval,1);
% %                     else
% %                         ruleval = rulevals;
% %                     end
% %                     nrule = SplitRule(X,colind,ruleval);
% %                 end
%             elseif strcmp(obj.Xclass(colind),'cell')
%                 S = unique(X(node.Xind,colind)); % unique groups
%                 ns = size(S,1);
%                 if ns == 1
%                     error('Need more than one group to split on')
%                 end    
%                 % Randomly Select a Subset
%                 nrule = SplitRule(X,colind,randsub(S));
%             else
%                 error('Unexpected class encountered.')
%             end
%         end
        
        function out = loglikefunc(obj,thetree,y)
            % LGP Implementation
%             m = thetree.m;
%             % Set optim options
%             opt = optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');  
%             % Find data which are in the node
%             ypart = y(obj.Xind);
%     
%             % Choose a subset of bins (to avoid estimating densities in regions
%             %   that don't have any data).
%             ymin = min([min(ypart) mean(ypart) - 3*std(ypart)]);
%             ymax = max([max(ypart) mean(ypart) + 3*std(ypart)]);
%             ylb = find(thetree.ygrid > ymin, 1 ) - 1;
%             if ylb < 1
%                 ylb = 1;
%             end
%             yup = find(thetree.ygrid < ymax, 1, 'last' ) + 1;
%             if yup > m
%                 yup = m;
%             end
%             ygridpart = thetree.ygrid(ylb:yup);
% 
%             % Calculate the number in each bin
%             nypart = hist(ypart,ygridpart);
%             nypart = nypart';
%             % Optimization of hyperparameters...
%             % Scale ygridpart so we can always use the same prior
%             % (same strategy is employed in lgpdens just before using gpsmooth
%             % function
%             %Xscaled = ((Xpart - mx)./sx)';
%             ygridpartn = ((ygridpart - mean(ygridpart))./std(ygridpart))';
%             % Best guess of lengthscale (From gpsmooth function in lgpdens from
%             %    riihimaki's excellent MATLAB code)
%             h=max(diff(ygridpartn(1:2,end)).^2,1/sum(nypart).^(1/5)/2);
%             gp = thetree.GP;
%             gpcf1 = gp.cf{1};
%             gpcf1 = gpcf_sexp(gpcf1, 'lengthScale', h*repmat(2,[1 size(ygridpartn,2)]));
%             gp = gp_set(gp,'cf',gpcf1);
%             gp = gp_optim(gp,ygridpartn,nypart,'opt',opt, 'optimf', @fminlbfgs);
%             % Change GP structure to drop the priors for the 'l' and 'sigma2'
%             %   The LGP paper just optimizes the hyperparameters and then
%             %   just uses those values.  The priors for the hyperparameters are
%             %   simply to aid in finding the smoothness parameters.  The marginal 
%             %   distribution in the paper is given the MAP hyperparameter values.
%             gp.cf{1} = gpcf_sexp(gp.cf{1},'magnSigma2_prior',[],'lengthScale_prior',[]);
%             % Calculate Marginal likelihood and add it to the rest.
%             % X or Xscaled for this part? I think Xscaled is fine and it is better
%             % numerically.  Also, we use Xscaled above, so it must be consistent
%             % with previous estimation.
%             % ll = ll - gpla_e([],gp,'x',X','y',nypart');
%             out = obj;
%             out.Llike = - gpla_e([],gp,'x',ygridpartn,'y',nypart);
            
            % Simple CART Implementation
            xind = obj.Xind;
            n = length(xind);
            alpha= 10;
            beta = var(y);
            yi = y(xind);
            mubar = mean(y);
            gma = sum(yi.^2) + mubar^2 - (sum(yi) + mubar)^2/(n + 1);
            llike = -n/2*log(2*pi) - .5*log(n + 1) + n/2*log(2) + ...
                gammaln((n+2)/2) - gammaln(alpha/2) + ...
                alpha/2 * log(beta) - (alpha + n)/2 * log(beta + gma);
            out = obj;
            out.Llike = llike;
        end
        
        function out = getsplits(obj,X,leafmin)
            out = obj;
            p = size(X,2);
            Xsub = X(obj.Xind,:);
            nsub = length(obj.Xind);
            % Determine which variables are available for a split
            %   and how many split points are available.
            nSplit = zeros(p,1);
            splitvals = cell(p,1);
            for ii=1:p
                xsub = Xsub{:,ii};
                nsplits = 0;
                svals = [];
                if isa(xsub,'numeric')
                    if length(unique(xsub)) > 1
                        thetab = tabulate(xsub);
                        % Remove any zero counts (can happen with only
                        %     positive integer x-variables
                        thetab = thetab(thetab(:,2) > 0,:);
                        
                        % Use the cumulative sums to figure out what split
                        %   variables are possible
                        csum = cumsum(thetab(:,2));
                        %rcsum = cumsum(thetab(:,2),'reverse');
                        % Index on the table of which values you can split
                        % on.
                        splitind = (csum >= leafmin) & ((nsub - csum) >= leafmin);
                        nsplits = sum(splitind);
                        svals = thetab(splitind,1);
                        if min(size(svals)) == 0
                            svals = [];
                        end
                    end
                elseif isa(xsub,'cell')
                    thetab = tabulate(xsub);
                    if length(unique(xsub)) > 1
                        % thetab = tabulate(xsub);
                        % group = thetab(:,1);
                        ngroup = size(thetab,1);
                        gdiv2 = floor(ngroup/2); % ngroup/2 and rounded down
%                         if all(cell2mat(thetab(:,2)) > leafmin) % Easy Case
%                             ncombn = zeros(gdiv2,1);
%                             for jj = 1:gdiv2
%                                 ncombn(jj) = nchoosek(ngroup,jj);
%                             end
%                             nsplits = sum(ncombn);
%                         else
                        % More comptuationally demanding case
                        for jj = 1:gdiv2
                            cmat = combnk(1:ngroup,jj);
%                             if jj == gdiv2 && mod(ngroup,2) == 0 % subset if necessary to avoid double counts
%                                 % ncmat = size(cmat,1);
%                                 I_cmat = cmat(:,1) == 1;
%                                 cmat = cmat(I_cmat,:);
%                             end
                            for kk = 1:size(cmat,1)
                                ngroup1 = sum(cell2mat(thetab(cmat(kk,:),2)));
                                ngroup2 = nsub - ngroup1;
                                if ngroup1 > leafmin && ngroup2 > leafmin
                                    nsplits = nsplits + 1;
                                    svals{nsplits} = thetab(cmat(kk,:),1);
                                end
                            end
                            if jj == gdiv2 && mod(ngroup,2) == 0
                                % In this case, we have "duplicate" rules
                                ncmat = size(cmat,1);
                                nsplits = nsplits - ncmat/2;
                            end
                        end
                        %end
                    elseif length(unique(xsub)) == 1 %&& all(cell2mat(thetab(:,2)) >= leafmin)
                        % nsplits = 1;
                        % svals{1} = thetab(1,1);
                    %elseif length(unique(xsub)) == 1 && ~all(cell2mat(thetab(:,2)) > leafmin)
                        % Use default values of nsplits an&& all(cell2mat(thetab(:,2)) > leafmin)d svals
                    else
                        c1 = length(unique(xsub))
                        c2 = all(cell2mat(thetab(:,2)) > leafmin)
                        thetab
                        % orig = obj
                        % newone = out
                        error('Unexpected case.')
                    end
                else
                    error('Data type received was not expected.')
                end
                
                % TODO: remove this case
                if nsplits == 1 && isempty(svals)
                    error('problem here')
                end
                
                nSplit(ii) = nsplits;
                splitvals{ii} = svals;
            end
            out.nSplits = nSplit;
            out.Splitvals = splitvals;
            out.Updatesplits = 0;
        end
        
        
    end
    methods(Static)
        % Create the prediction matrix
        % Returns a px1 vector where the ith element indicates the number
        %   of possible splits for the ith predictor.
%         function out = createpredmatrix(X,Xind,Leafmin)
%             p = size(X,2);
%             Xsub = X(Xind,:);
%             nsub = length(Xind);
%             % Determine which variables are available for a split
%             %   and how many split points are available.
%             out = zeros(p,1);
%             for ii=1:p
%                 xsub = Xsub{:,ii};
%                 nsplits = 0;
%                 if isa(xsub,'numeric')
%                     if length(unique(xsub)) > 1
%                         thetab = tabulate(xsub);
%                         % Remove any zero counts (can happen with only
%                         %     positive integer x-variables
%                         thetab = thetab(thetab(:,2) > 0,:);
%                         
%                         % Use the cumulative sums to figure out what split
%                         %   variables are possible
%                         csum = cumsum(thetab(:,2));
%                         %rcsum = cumsum(thetab(:,2),'reverse');
%                         % Index on the table of which values you can split
%                         % on.
%                         splitind = (csum >= Leafmin) & ((nsub - csum) >= Leafmin);
%                         nsplits = sum(splitind);
%                     end
%                 elseif isa(xsub,'cell')
%                     if length(unique(xsub)) > 2
%                         thetab = tabulate(xsub);
%                         % group = thetab(:,1);
%                         ngroup = size(thetab,1);
%                         gdiv2 = floor(ngroup/2); % ngroup/2 and rounded down
%                         if all(cell2mat(thetab(:,2)) > Leafmin) % Easy Case
%                             ncombn = zeros(gdiv2,1);
%                             for jj = 1:gdiv2
%                                 ncombn(jj) = nchoosek(ngroup,jj);
%                             end
%                             nsplits = sum(ncombn);
%                         else
%                             % More comptuationally demanding case
%                             for jj = 1:gdiv2
%                                 cmat = combnk(1:ngroup,jj);
%                                 for kk = 1:size(cmat,1)
%                                     ngroup1 = sum(cell2mat(thetab(cmat(kk,:),2)));
%                                     ngroup2 = nsub - ngroup1;
%                                     if ngroup1 > Leafmin && ngroup2 > Leafmin
%                                         nsplits = nsplits + 1;
%                                     end
%                                 end
%                             end
%                         end
%                     elseif length(unique(xsub)) == 2 && all(cell2mat(thetab(:,2)) > Leafmin)
%                         nsplits = 1;
%                     end
%                 else
%                     error('Data type received was not expected.')
%                 end
%                 out(ii) = nsplits;
%             end
%         end
    end
end