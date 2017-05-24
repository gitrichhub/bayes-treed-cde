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


function [output] = TreeMCMC_Multiset(y,X,nmcmc,burn,K,leafmin,gamma,beta,p)
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
    
    % Initialize root tree
    T = Tree(y,X,leafmin,gamma,beta);
    %startprior = prior_eval(T,X);
    % Make copies of the root tree for the multiset of size K
    Tmultiset = cell(K,1);
    %Tpriors = zeros(K,1);
    for ii = 1:K
        Tmultiset{ii} = T;
        %Tpriors(ii) = startprior;
    end
    % Remove after debugging
    T = [];
    %startprior = [];
    
    % Probability of proposing steps
    p_g_orig = .25; % grow
    p_p_orig = .25; % prune
    p_c_orig = .25; % change
    p_s_orig = .25; % swap
    allprobs = [p_g_orig,p_p_orig,p_c_orig,p_s_orig];
    allprobs = allprobs./sum(allprobs);
    
    n_totals = zeros(4,1); % total tries for grow, prune, change, & swap steps, respectively.

    naccept = 0;
%     n_g_accept = 0;
%     n_p_accept = 0;
%     n_c_accept = 0;
%     n_s_accept = 0;
    
    
    % Initial Values
    post = getpost(Tmultiset,X);
    post = ones(K,1)*post;
    
    
    % Posterior trees multisets
    TREEMULTISETS = cell(nmcmc,1);
    treesize = zeros(nmcmc,K);
    LLIKE = zeros(nmcmc,K);
    tsizes = ones(1,K);
    moveaccept = zeros(nmcmc,7);
    for ii=1:(burn + nmcmc)
        accept = 0;
        % Randomly choose a member of the multiset to update
        kind = randsample(K,1);
        [Tstar,~,n_totals,prop_ratio,tsize,r,~] = proposeTree(Tmultiset{kind},[],y,X,allprobs,n_totals,p,1);
        %problem = getsplitcheck(Tstar,X);
        %prop_ratio = 0; % THIS IS ONLY EXPERIMENTAL! REMOVE THIS LINE!
        Tmultisetstar = Tmultiset;
        Tmultisetstar{kind} = Tstar;
        % Obtain multiplicity
        m1 = mfunc(kind,Tmultiset);
        m2 = mfunc(kind,Tmultisetstar);
        [poststar,llikestar,~] = getpost(Tmultisetstar,X);
        postratio = (poststar - post(kind));
        lr = postratio + log(m2/m1) + prop_ratio;
        
        
        if lr > log(rand)
            accept = 1;
            Tmultiset = Tmultisetstar;
            post(kind) = poststar;
            llikes = llikestar;
            naccept = naccept + 1;
            tsizes(kind) = tsize;
%             if r == 1
%                 n_g_accept = n_g_accept + 1;
%                 tsize = tsize + 1;
%             elseif r == 2
%                 n_p_accept = n_p_accept + 1;
%                 tsize = tsize - 1;
%             elseif r == 3
%                 n_c_accept = n_c_accept + 1;
%             else
%                 n_s_accept = n_s_accept + 1;
%             end
        end
        
        if mod(ii,10) == 0
            disp(['i = ',num2str(ii),', llike = ',num2str(llikes'),...
                ...%', lliketotal = ',num2str(sum(llikes)),...
                ', accept = ',num2str(naccept/ii),...
                ', Size = ',num2str(tsizes)]);
            % Treeplot(T);
        end
        
        
        if ii > burn
            TREEMULTISETS{ii - burn} = Tmultiset;
            treesize(ii - burn,:) = tsizes;
            LLIKE(ii - burn,:) = llikes;    
            moveaccept(ii-burn,:) = [kind,r,accept,postratio,log(m2/m1),prop_ratio,lr];
        end
    end
    perc_accept = naccept/(nmcmc + burn);
%     move_accepts = [n_g_accept/n_g_total,...
%         n_p_accept/n_p_total,...
%         n_c_accept/n_c_total,...
%         n_s_accept/n_s_total];
    moveacceptfinal = array2table(moveaccept,'VariableNames',{'K_ind','move','accepted','loglikeratio','logmultiplicityratio',...
        'prop_ratio','lr'});
    output = struct('Treemultisets',{TREEMULTISETS},'llike',LLIKE,'acceptance',perc_accept,...
        'treesize',treesize,'moveaccept',moveacceptfinal);%,'move_accepts',move_accepts);
end


