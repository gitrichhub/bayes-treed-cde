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


function [output] = TreeMCMCtemp(y,X,nmcmc,burn,leafmin,gamma,beta,p,temp)
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
    T = Tree(y,X,leafmin,gamma,beta,temp);
    
    % Probability of proposing steps
    p_g_orig = .25; % grow
    p_p_orig = .25; % prune
    p_c_orig = .25; % change
    p_s_orig = .25; % swap
    allprobs = [p_g_orig,p_p_orig,p_c_orig,p_s_orig];

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
        [Tstar,~,r,lr] = proposeTree(T,y,X,allprobs,p,temp,0);
        if lr > log(rand)
            T = Tstar;
            naccept = naccept + 1;
            if r == 1
                n_g_accept = n_g_accept + 1;
                n_g_total = n_g_total + 1;
                tsize = tsize + 1;
            elseif r == 2
                n_p_accept = n_p_accept + 1;
                n_p_total = n_p_total + 1;
                tsize = tsize - 1;
            elseif r == 3
                n_c_accept = n_c_accept + 1;
                n_c_total = n_c_total + 1;
            else
                n_s_accept = n_s_accept + 1;
                n_s_total = n_s_total + 1;
            end
        else
            if r == 1
                n_g_total = n_g_total + 1;
            elseif r == 2
                n_p_total = n_p_total + 1;
            elseif r == 3
                n_c_total = n_c_total + 1;
            else
                n_s_total = n_s_total + 1;
            end
        end
        
        if mod(ii,10) == 0
            disp(['i = ',num2str(ii),', llike = ',num2str(T.Lliketree),...
                ', accept = ',num2str(naccept/ii),...
                ', Size = ',num2str(tsize)]);
        end
        % Record Values       
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

    