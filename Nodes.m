classdef Nodes
    properties
        Id
        Parent
        Lchild
        Rchild
        Rule
        Xind % index of which data values pass through if a terminal node
        Llike % Log-Likelihood contribution of data in the partition...
        Updatellike % 0 if no llike computation needed, 1 if it is needed
    end
    methods
        % Constructor
        function obj = Nodes(Id,Parent,Lchild,Rchild,Rule,Xind) % Add in an Xind later?
            if isa(Id,'numeric')
                obj.Id = Id;
            else
                error('Id must be numeric')
            end
            obj.Parent = Parent;
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
            if isa(Rule,'SplitRule')
                obj.Rule = Rule;
            elseif isempty(Rule)
            else
                error('Rule must be of class "SplitRule" or empty')
            end
            
            if isa(Xind,'numeric')
                obj.Xind = Xind;
            elseif isempty(Xind)
            else
                error('Xind must be numeric or empty')
            end
            % Llike needs to be computed.
            obj.Updatellike = 1;
        end
        
        % Add/Change Rule
        function out = newrule(obj,therule)
            if isa(therule,'SplitRule') || isempty(therule)
                obj.Rule = therule;
                out = obj;
            else
                error('rule must be of class "SplitRule"');
            end
        end
        
        % select rule for factors
        function out = factrule(obj)
            out.Xind
        end
        
       
        
        
        
        
        
        
        % Draw a new rule from the prior
        function [nrule, colind] = drawrule(node,obj,X,colind) % obj is a Tree class object
            % Choosing colind = 0 will randomly choose a variable to split
            % Randomly choose a column to split on
            if colind < 1
                % DO SOME CHECKING BEFORE CHOOSING THE COLUMN:
                %   check to see if there are at least two categories
                %   check to make sure we have enough data to justify the
                %   split.
                % Randomly choose a column to split on (if not specified)
                candcols = [];
                for ii = 1:size(X,2)
                    %if strcmp(obj.Xclass(ii),'double')
                        nu = max(size(unique(X(node.Xind,ii))));
                        if nu > 1
                            candcols = [candcols, ii];
                        end
                    %elseif strcmp(obj.Xclass(ii),'cell')
                end
                if ~isempty(candcols)
                    colind = randsample(candcols,1);
                else
                    nrule = [];
                    colind = 0;
                    return
                end
            end
            if strcmp(obj.Xclass(colind),'double')
                rulevals = randsample(table2array(X(node.Xind,colind)),length(node.Xind));
                
                % OLD VERSION...
                % Don't allow the rule to be the maximum (gives an empty
                % node)
                mval = max(rulevals);
                for ii = 1:length(rulevals)
                    if rulevals(ii) < mval
                        nrule = SplitRule(X,colind,rulevals(ii));
                        break
                    end
                end
                

                % New Version
                % Subset rules so that you have at least leafmin obs at
                %   each terminal node...
%                 rulevals = sort(table2array(X(node.Xind,colind)));
%                 len = length(rulevals);
%                 leafmin = obj.Leafmin;
%                 if len >= 2*leafmin
%                     rulevals = rulevals([leafmin,len-leafmin]);
%                     if length(rulevals) > 1
%                         ruleval = randsample(ruleval,1);
%                     else
%                         ruleval = rulevals;
%                     end
%                     nrule = SplitRule(X,colind,ruleval);
%                 end
            elseif strcmp(obj.Xclass(colind),'cell')
                S = unique(X(node.Xind,colind)); % unique groups
                ns = size(S,1);
                if ns == 1
                    error('Need more than one group to split on')
                end    
                % Randomly Select a Subset
                nrule = SplitRule(X,colind,randsub(S));
            else
                error('Unexpected class encountered.')
            end
        end
        
        function out = loglikefunc(obj,thetree,y)
            m = thetree.m;
            % Set optim options
            opt = optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');  
            % Find data which are in the node
            ypart = y(obj.Xind);
    
            % Choose a subset of bins (to avoid estimating densities in regions
            %   that don't have any data).
            ymin = min([min(ypart) mean(ypart) - 3*std(ypart)]);
            ymax = max([max(ypart) mean(ypart) + 3*std(ypart)]);
            ylb = find(thetree.ygrid > ymin, 1 ) - 1;
            if ylb < 1
                ylb = 1;
            end
            yup = find(thetree.ygrid < ymax, 1, 'last' ) + 1;
            if yup > m
                yup = m;
            end
            ygridpart = thetree.ygrid(ylb:yup);

            % Calculate the number in each bin
            nypart = hist(ypart,ygridpart);
            nypart = nypart';
            % Optimization of hyperparameters...
            % Scale ygridpart so we can always use the same prior
            % (same strategy is employed in lgpdens just before using gpsmooth
            % function
            %Xscaled = ((Xpart - mx)./sx)';
            ygridpartn = ((ygridpart - mean(ygridpart))./std(ygridpart))';
            % Best guess of lengthscale (From gpsmooth function in lgpdens from
            %    riihimaki's excellent MATLAB code)
            h=max(diff(ygridpartn(1:2,end)).^2,1/sum(nypart).^(1/5)/2);
            gp = thetree.GP;
            gpcf1 = gp.cf{1};
            gpcf1 = gpcf_sexp(gpcf1, 'lengthScale', h*repmat(2,[1 size(ygridpartn,2)]));
            gp = gp_set(gp,'cf',gpcf1);
            gp = gp_optim(gp,ygridpartn,nypart,'opt',opt, 'optimf', @fminlbfgs);
            % Change GP structure to drop the priors for the 'l' and 'sigma2'
            %   The LGP paper just optimizes the hyperparameters and then
            %   just uses those values.  The priors for the hyperparameters are
            %   simply to aid in finding the smoothness parameters.  The marginal 
            %   distribution in the paper is given the MAP hyperparameter values.
            gp.cf{1} = gpcf_sexp(gp.cf{1},'magnSigma2_prior',[],'lengthScale_prior',[]);
            % Calculate Marginal likelihood and add it to the rest.
            % X or Xscaled for this part? I think Xscaled is fine and it is better
            % numerically.  Also, we use Xscaled above, so it must be consistent
            % with previous estimation.
            % ll = ll - gpla_e([],gp,'x',X','y',nypart');
            out = obj;
            out.Llike = - gpla_e([],gp,'x',ygridpartn,'y',nypart);
        end
    end
    
end