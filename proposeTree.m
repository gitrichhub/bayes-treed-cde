function [Tstar,prop_ratio,r,lr] = proposeTree(T,y,X,allprobs,p,temp,mset)
    % T: current tre
    % y: dependent variable
    % X: Data
    Tprior = T.Prior;
    
    % Initial values
    lr = [];
    %n_totals = ntotals;

    % See what moves are possible
    if length(T.Allnodes) >= 5
        [~,swappossible] = swap(T,y,X,[],0);
    else
        swappossible = [];
    end
    [nbirths,T] = nbirthnodes(T,X);
    [p_g,p_p,p_c,p_s] = propprobs(T,allprobs,nbirths,swappossible);
    r = randsample([1, 2, 3, 4],1,true,[p_g,p_p,p_c,p_s]);

    if r == 1 % grow
        [Tstar,birthindex] = birth(T,y,X);
        % kstar = length(termnodes(T)); % number of birthable nodes
        kstar = nbirthnodes(T,X); % Number of nodes that can grow
        k_d = length(terminalparents(Tstar)); % number of possible deaths in proposed model
        birthvarind = Tstar.Allnodes{birthindex}.Rule{1};
        % Number of possible variables to split on
        N_v = sum(Tstar.Allnodes{birthindex}.nSplits > 0);
        % Number of possible splits for this variable
        N_b = Tstar.Allnodes{birthindex}.nSplits(birthvarind);

        % Reversibility
        if length(Tstar.Allnodes) >= 5
            [~,swappossible] = swap(Tstar,y,X,[],0);
        else
            swappossible = [];
        end
        [nbirths,Tstar] = nbirthnodes(Tstar,X);
        [~,p_p_star,~,~] = propprobs(Tstar,allprobs,nbirths,swappossible);

        prop_ratio = log(p_p_star/k_d) - log(p_g/(kstar*N_v*N_b));
        [Tstarprior,Tstar] = prior_eval(Tstar,X);
        
        if ~mset
            lr = temp*(Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - Tprior) + ...
                prop_ratio;
        end
        %n_g_total = n_g_total + 1;
        % n_totals(1) = ntotals(1) + 1;
    elseif r == 2 % prune
        [Tstar,pind] = prune(T,y,X);
        %kstar = length(termnodes(Tstar));
        k_d = length(terminalparents(T));
        prunevarind = T.Allnodes{pind}.Rule{1};
        N_b = T.Allnodes{pind}.nSplits(prunevarind);
        N_v = sum(T.Allnodes{pind}.nSplits > 0);
        % Reversibility
        if length(Tstar.Allnodes) >= 5
            [~,swappossible] = swap(Tstar,y,X,[],0);
        else
            swappossible = [];
        end
        [nbirths,Tstar] = nbirthnodes(Tstar,X);
        kstar = nbirths;
        [p_g_star,~,~,~] = propprobs(Tstar,allprobs,nbirths,swappossible);

        prop_ratio = log(p_g_star/(kstar*N_v*N_b)) - log(p_p/k_d);
        [Tstarprior,Tstar] = prior_eval(Tstar,X);
        if ~mset
            lr = temp*(Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - Tprior) + ...
                prop_ratio;
        end
        %n_p_total = n_p_total + 1;
        % n_totals(2) = ntotals(2) + 1;
    elseif r == 3; % change
        [Tstar,priordraw,startcont,endcont,nchange,nchange2] = change(T,y,X,p);
        % Tstar = change(T,y,X);
        [Tstarprior,Tstar] = prior_eval(Tstar,X);
        %[nT,T] = nchangenodes(T,X);
        %[nTstar,Tstar] = nchangenodes(Tstar,X);

        % Reversibility
        if priordraw
            if startcont
                if endcont            
                    prop_ratio = 0;
                else
                    prop_ratio = log(1/(1-p));
                end
            else % start with categorical variable
                if endcont
                    prop_ratio = log(1-p);
                else
                    prop_ratio = 0;
                end
            end
        else % not a prior draw
            if startcont
                if endcont
                    if nchange > 0 && nchange2 > 0
                        prop_ratio = log(nchange/nchange2);
                    elseif nchange == 0 && nchange2 == 0
                        prop_ratio = 0;
                    else
                        error('Should not occur')
                    end
                else
                    error('This should not happen.')
                end
            else
                error('This should never happen.')
            end
        end

        % prop_ratio = log(nT/nTstar);
        if ~mset
            lr = temp*(Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - Tprior) + ...
                prop_ratio;
        end
        %n_c_total = n_c_total + 1;
        %n_totals(3) = ntotals(3) + 1;
    elseif r == 4; % swap
        Tstar = swap(T,y,X,[],1);
        [Tstarprior,Tstar] = prior_eval(Tstar,X);
        nT = nswaps(T,y,X);
        nTstar = nswaps(Tstar,y,X);
        prop_ratio = log(nT/nTstar);
        if ~mset
            lr = temp*(Tstar.Lliketree - T.Lliketree) + ...
                Tstarprior - Tprior + prop_ratio;
        end
        %n_s_total = n_s_total + 1;
        %n_totals(4) = ntotals(4) + 1;
    end
    
    if ~isfinite(lr)
        error('Non-finite likelihood ratio')
    end
    
    %tsize = ceil(length(Tstar.Allnodes)/2);
end
    
    