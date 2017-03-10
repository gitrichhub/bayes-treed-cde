function [TREES,perc_accept] = TreeMCMC(y,X,nmcmc,burn,leafmin,gamma,beta)
    % Check for format of input variables
    if ~isnumeric(y)
        error('y must be numeric.')
    end
    if ~istable(X)
        error('X must be a table.')
    end
    if length(nmcmc) ~= 1 
        error('input arg "nmcmc" must be a scalar integer.')
    end
    if length(burn) ~= 1
        error('input arg "burn" must be a scalar integer.')
    end
    
    
    % Initialize root tree
    T = Tree(y,X,leafmin,gamma,beta);
    
    % Probability of birth and deaths
    p_b_orig = .5;
    p_d_orig = .5;
    
    % 
    naccept = 0;
    
    % Posterior trees
    TREES = cell(nmcmc,1);
    for ii=1:(burn + nmcmc)
        % Simple Grow or Prune functionality  

        % Choose grow or prune step
        if rand < p_b_orig
            r = 1;
        else
            r = 2;
        end
        % Adjustment for boundary case
        p_d = p_d_orig;
        if length(T.Allnodes) == 1
            r = 1;
            p_b = 1;
        else
            p_b = p_b_orig;
        end

        if r == 1 % grow
            [Tstar,birthindex] = birth(T,y,X);
            kstar = length(termnodes(T));
            k_d = length(terminalparents(Tstar)); % number of possible deaths in proposed model
            birthvarind = Tstar.Allnodes{birthindex}.Rule{1};
            % Number of possible splits for this variable
            N_b = Tstar.Allnodes{birthindex}.nSplits(birthvarind);
            % For simplicity, we assume all nodes are capable of
            %   growing.  However, this will need to be adjusted later
            %   when we have nodes which cannot be split further.
            prop_ratio = log(p_d/k_d) - log(p_b/(kstar*N_b));
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (prior_eval(Tstar) - prior_eval(T)) + ...
                prop_ratio;
        elseif r == 2 % prune
            [Tstar,pind] = prune(T,y);
            kstar = length(termnodes(Tstar));
            k_d = length(terminalparents(T));
            prunevarind = T.Allnodes{pind}.Rule{1};
            N_b = T.Allnodes{pind}.nSplits(prunevarind);
            prop_ratio = log(p_b/(kstar*N_b)) - log(p_d/k_d);
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (prior_eval(Tstar) - prior_eval(T)) + ...
                prop_ratio;
        end

        if lr > log(rand)
            T = Tstar;
            naccept = naccept + 1;
        end
        
        if mod(ii,10) == 1
            disp(['i = ',num2str(ii),', llike = ',num2str(T.Lliketree)])
            Treeplot(T);
        end
        
        if ii > burn
            TREES{ii - burn} = T;
        end
    end
    perc_accept = naccept/(nmcmc + burn);
end


