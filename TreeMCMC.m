function [output] = TreeMCMC(y,X,nmcmc,burn,leafmin,gamma,beta)
    % TODO: 
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
    
    % Probability of proposing steps
    p_g_orig = .25; % grow
    p_p_orig = .25; % prune
    p_c_orig = .25; % change
    p_s_orig = .25; % swap
    allprobs = [p_g_orig,p_p_orig,p_c_orig,p_s_orig];
%     cum_prob = cumsum([p_g_orig,p_p_orig,p_c_orig,p_s_orig]);
%     cum_prob2 = cumsum([p_g_orig,p_p_orig,p_c_orig])./ ...
%         (p_g_orig+p_p_orig+p_c_orig);
%     tot_prob = p_g_orig + p_p_orig + p_c_orig;
    
    % 
    naccept = 0;
    n_g_accept = 0;
    n_p_accept = 0;
    n_c_accept = 0;
    n_s_accept = 0;
    n_g_total = 0;
    n_p_total = 0;
    n_c_total = 0;
    n_s_total = 0;
    
    % Posterior trees
    TREES = cell(nmcmc,1);
    treesize = zeros(nmcmc,1);
    LLIKE = zeros(nmcmc,1);
    tsize = 1;
    for ii=1:(burn + nmcmc)
        % See what moves are possible
        if length(T.Allnodes) >= 5
            [~,swappossible] = swap(T,y,X,[]);
        else
            swappossible = [];
        end
        [nbirths,T] = nbirthnodes(T,X);
        [p_g,p_p,p_c,p_s] = propprobs(T,allprobs,nbirths,swappossible);
        r = randsample([1, 2, 3, 4],1,true,[p_g,p_p,p_c,p_s]);
        
%         
%         
%         if nbirths > 0 && ~isempty(swappossible)
%             rval = rand;
%             % Simple Grow or Prune functionality  
%             if length(T.Allnodes) >= 5 % all steps are possible (note that 
%                 %                        a swap step requires at least 5 nodes)
%                 % Choose grow, prune, change, or swap steps
%                 if rval < cum_prob(1) % grow
%                     r = 1;
%                 elseif rval < cum_prob(2) % prune
%                     r = 2; 
%                 elseif rval < cum_prob(2) % change
%                     r = 3;
%                 else % swap
%                     r = 4;
%                 end
%             elseif length(T.Allnodes) > 1
%                 if rval < cum_prob2(1) % grow
%                     r = 1;
%                 elseif rval < cum_prob2(2) % prune
%                     r = 2;
%                 else 
%                     r = 3; % change
%                 end
%             else
%                 r = 1;
%             end
%             % Adjustment for MH for boundary cases
%             p_g = p_g_orig;
%             p_p = p_p_orig;
%             %p_c = p_c_orig;
%             %p_s = p_s_orig;
%             if length(T.Allnodes) == 1
%                 p_g = 1;
%                 p_p = p_p_orig/tot_prob;
%             elseif length(T.Allnodes) == 4
%                 if r == 1
%                     p_g = p_g_orig/tot_prob;
%                     % p_p = p_p_orig; % already defined above
%                 elseif r == 2
%                     p_g = p_g_orig/tot_prob;
%                     p_p = p_p_orig/tot_prob;
%                     %p_c = p_c_orig/tot_prob;
%                 end
%             elseif length(T.Allnodes) == 5
%                 if r == 1 % no change necessary
%                 elseif r == 2
%                     % p_p = p_p_orig;
%                     p_g = p_g_orig/tot_prob;
%                 end
%             end
%         else
%             moves = [2,3]; % possible move steps
%             if nbirths > 0
%                 moves = [1,moves];
%                 p_g
%             end
%             if ~isempty(swappossible)
%                 moves = [moves,4];
%             end
%             probs = allprobs(moves);
%             probs = probs./sum(probs);
%             r = randsample(moves,1,true,probs);
            
                

        if r == 1 % grow
            [Tstar,birthindex] = birth(T,y,X);
            % kstar = length(termnodes(T)); % number of birthable nodes
            kstar = nbirthnodes(T,X); % Number of nodes that can grow
            k_d = length(terminalparents(Tstar)); % number of possible deaths in proposed model
            birthvarind = Tstar.Allnodes{birthindex}.Rule{1};
            % Number of possible splits for this variable
            N_b = Tstar.Allnodes{birthindex}.nSplits(birthvarind);
            
            % Reversibility
            if length(Tstar.Allnodes) >= 5
                [~,swappossible] = swap(Tstar,y,X,[]);
            else
                swappossible = [];
            end
            [nbirths,Tstar] = nbirthnodes(Tstar,X);
            [~,p_p_star,~,~] = propprobs(Tstar,allprobs,nbirths,swappossible);
            
            prop_ratio = log(p_p_star/k_d) - log(p_g/(kstar*N_b));
            [Tstarprior,Tstar] = prior_eval(Tstar,X);
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - prior_eval(T,X)) + ...
                prop_ratio;
            n_g_total = n_g_total + 1;
        elseif r == 2 % prune
            [Tstar,pind] = prune(T,y);
            %kstar = length(termnodes(Tstar));
            k_d = length(terminalparents(T));
            prunevarind = T.Allnodes{pind}.Rule{1};
            N_b = T.Allnodes{pind}.nSplits(prunevarind);
            
            % Reversibility
            if length(Tstar.Allnodes) >= 5
                [~,swappossible] = swap(Tstar,y,X,[]);
            else
                swappossible = [];
            end
            [nbirths,Tstar] = nbirthnodes(Tstar,X);
            kstar = nbirths;
            [p_g_star,~,~,~] = propprobs(Tstar,allprobs,nbirths,swappossible);
            
            prop_ratio = log(p_g_star/(kstar*N_b)) - log(p_p/k_d);
            [Tstarprior,Tstar] = prior_eval(Tstar,X);
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - prior_eval(T,X)) + ...
                prop_ratio;
            n_p_total = n_p_total + 1;
        elseif r == 3; % change
            Tstar = change(T,y,X);
            % Assuming the number of possible changes are the same for reversibility
            %   this is not necessarily true, but is difficult to compute.
            %   Perhaps will add this in future version
            [Tstarprior,Tstar] = prior_eval(Tstar,X);
            [nT,T] = nchangenodes(T,X);
            [nTstar,Tstar] = nchangenodes(Tstar,X);
            
             % Reversibility (nothing needed for change step)           
            
            prop_ratio = log(nT/nTstar);
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - prior_eval(T,X)) + ...
                prop_ratio;
            n_c_total = n_c_total + 1;
        elseif r == 4; % swap
            Tstar = swap(T,y,X,[]);
            [Tstarprior,Tstar] = prior_eval(Tstar,X);
            nT = nswaps(T,y,X);
            nTstar = nswaps(Tstar,y,X);
            prop_ratio = log(nT/nTstar);
            lr = Tstarprior - prior_eval(T,X) + prop_ratio;
            n_s_total = n_s_total + 1;
        end

        if lr > log(rand)
            T = Tstar;
            naccept = naccept + 1;
            if r == 1
                n_g_accept = n_g_accept + 1;
                tsize = tsize + 1;
            elseif r == 2
                n_p_accept = n_p_accept + 1;
                tsize = tsize - 1;
            elseif r == 3
                n_c_accept = n_c_accept + 1;
            else
                n_s_accept = n_s_accept + 1;
            end
        end
        
        if mod(ii,10) == 0
            disp(['i = ',num2str(ii),', llike = ',num2str(T.Lliketree),...
                ', accept = ',num2str(naccept/ii),...
                ', Size = ',num2str(tsize)]);
            % Treeplot(T);
        end
        
        
        if ii > burn
            TREES{ii - burn} = T;
            treesize(ii - burn) = tsize;
            LLIKE(ii - burn) = T.Lliketree;
            
        end
    end
    perc_accept = naccept/(nmcmc + burn);
    move_accepts = [n_g_accept/n_g_total,...
        n_p_accept/n_p_total,...
        n_c_accept/n_c_total,...
        n_s_accept/n_s_total];
        
    output = struct('Trees',{TREES},'llike',LLIKE,'acceptance',perc_accept,...
        'treesize',treesize,'move_accepts',move_accepts);
end

% T: an object of class T
% probs: a vector of the original proposal probabilities of the steps
% nbirths: number of births the Tree has available
% swappossible: [] or 1 indicating no/yes to a swap step on tree T
function [p_g,p_p,p_c,p_s] = propprobs(T,probs,nbirths,swappossible)
    if length(T.Allnodes) >= 5
        if nbirths > 0 % birth step possible
            if ~isempty(swappossible) % If swap is possible
                p_g = probs(1);
                p_p = probs(2);
                p_c = probs(3);
                p_s = probs(4);
            else % swap not possible
                prob_t = probs(1:3)./sum(probs(1:3));
                p_g = prob_t(1);
                p_p = prob_t(2);
                p_c = prob_t(3);
                p_s = 0;
            end
        else % If birth step is not possible
            if ~isempty(swappossible)  % swap possible
                prob_t = probs(2:4)./sum(probs(2:4));
                p_g = 0;
                p_p = prob_t(1);
                p_c = prob_t(2);
                p_s = prob_t(3);
            else % swap not possible
                prob_t = probs(2:3)./sum(probs(2:3));
                p_g = 0;
                p_p = prob_t(1);
                p_c = prob_t(2);
                p_s = 0;
            end
        end
    elseif length(T.Allnodes) > 1 % Only grow, death, and change steps available
        p_s = 0;
        if nbirths > 0 % grow step possible
            prob_t = probs(1:3)./sum(probs(1:3));
            p_g = prob_t(1);
            p_p = prob_t(2);
            p_c = prob_t(3);
        else % If grow step is not possible
            prob_t = probs(2:3)./sum(probs(2:3));
            p_g = 0;
            p_p = prob_t(1);
            p_c = prob_t(2);
        end
    else % root node
        p_g = 1;
        p_p = 0;
        p_c = 0;
        p_s = 0;
    end
end
    

    