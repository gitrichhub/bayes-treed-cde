%     This file is part of bayes-treed-cde.
% 
%     bayes-treed-cde is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     bayes-treed-cde is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with bayes-treed-cde.  If not, see <http://www.gnu.org/licenses/>.
%
%     Copyright 2016-2017, Richard Payne


function [output] = TreeMCMC(y,X,nmcmc,burn,leafmin,gamma,beta,p)
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
    
    if isempty(p)
        p = .75;
    end
    
    
%     for ii = 1:size(X,2)
%         if isa(X{:,ii},'cell')
%             break
%         end
%         if ii == size(X,2) % if all X's are continuous
%             p = 1;
%             disp('NOTE: p changed to 1 since all Xs are continuous');
%         end
%     end
            
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
    Tprior = prior_eval(T,X);
    for ii=1:(burn + nmcmc)
        % See what moves are possible
        if length(T.Allnodes) >= 5
            [~,swappossible] = swap(T,y,X,[],0);
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
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - Tprior) + ...
                prop_ratio;
            n_g_total = n_g_total + 1;
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
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - Tprior) + ...
                prop_ratio;
            n_p_total = n_p_total + 1;
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
                        prop_ratio = log(nchange/nchange2);
                    else
                        error('This should not happen.')
                    end
                else
                    error('This should never happen.')
                end
            end
            
            % prop_ratio = log(nT/nTstar);
            lr = (Tstar.Lliketree - T.Lliketree) + ...
                (Tstarprior - Tprior) + ...
                prop_ratio;
            n_c_total = n_c_total + 1;
        elseif r == 4; % swap
            Tstar = swap(T,y,X,[],1);
            [Tstarprior,Tstar] = prior_eval(Tstar,X);
            nT = nswaps(T,y,X);
            nTstar = nswaps(Tstar,y,X);
            prop_ratio = log(nT/nTstar);
            lr = Tstarprior - Tprior + prop_ratio;
            n_s_total = n_s_total + 1;
        end

        if lr > log(rand)
            T = Tstar;
            Tprior = Tstarprior;
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

    