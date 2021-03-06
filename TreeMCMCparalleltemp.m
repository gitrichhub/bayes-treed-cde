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

% This function performs MCMC using parallel tempering to search the
%   posterior of trees for conditional density estimation. This 
%   function writes the results to a file which can then be loaded into the
%   workspace.
%
% Required inputs:
% y: response variable, vector of a continuous response.
% X: table of covariates
%
% Optional Arguments:
% nmcmc: number of MCMC iterations, default 10,000
% burn: number of burn in iterations, default 1,000
% leafmin: the minimum number of observations required at a terminal node.
%          Default is 25.
% gamma: prior hyperparameter governing size/shape of tree.  See Chipman's
%        CART paper for details. Default is .95.
% beta: prior hyperparameter governing size/shape of tree.  See Chipman's
%        CART paper for details. Default is 1.
% k: The number of posterior trees to be returned.  If 0, all trees from
%    the MCMC algorithm are returned.  Otherwise, the k trees with the
%    highest marginal likelihoods are returned to save on memory. Default
%    is 0.
% p: the probability of performing an incremental change step on a
%    continuous variable rather than drawing from the prior. Default is
%    .75.
% parallelprofile: If 1, parallel profiling is performed to determine code
%                  efficiency.  Default is 0.
% hottemp: the smallest inverse temperature to be used.  Default is .1
% saveall: If 1, the results of all tempered chains are saved as output.
%          Default is 0.
% swapfreq: The tempered MCMC chains are swapped every 'swapfeq' iterations.
%           Default is 1. 
% seed: specify a seed for the MCMC run.  This creates independent random seeds
%       across each of the tempered chains.  Default is a random seed,
%       'shuffle'.
% suppress_errors_on_workers: If 1, suppresses warnings on workers. Default
%                             is 0. 
% filepath: where to place MCMC output. Default is './output/'

function TreeMCMCparalleltemp(y,X,varargin)
    % Parse function
    ip = inputParser;
    ip.FunctionName = 'TreeMCMCparalleltemp';
    % Required Inputs
    addRequired(ip,'y',@isnumeric);
    addRequired(ip,'X',@istable);
    % Optional Inputs
    addParameter(ip,'nmcmc',10000)
    addParameter(ip,'burn',1000)
    addParameter(ip,'leafmin',25)
    addParameter(ip,'gamma',.95)
    addParameter(ip,'beta',1)
    addParameter(ip,'k',0)
    addParameter(ip,'p',.75)
    addParameter(ip,'parallelprofile',0);
    addParameter(ip,'hottemp',.1)
    addParameter(ip,'saveall',0);
    addParameter(ip,'swapfreq',1);
    addParameter(ip,'seed','shuffle');
    addParameter(ip,'suppress_errors_on_workers',0);
    addParameter(ip,'filepath','./output/')
    
   
    parse(ip,y,X,varargin{:});
    y = ip.Results.y;
    X = ip.Results.X;
    nmcmc = ip.Results.nmcmc;
    burn = ip.Results.burn;
    leafmin = ip.Results.leafmin;
    gamma = ip.Results.gamma;
    beta = ip.Results.beta;
    k = ip.Results.k;
    p = ip.Results.p;
    parallelprofile = ip.Results.parallelprofile;
    hottemp = ip.Results.hottemp;
    saveall = ip.Results.saveall;
    seed = ip.Results.seed;
    swapfreq = ip.Results.swapfreq;
    suppress_errors_on_workers = ip.Results.suppress_errors_on_workers;
    filepath = ip.Results.filepath;
    
    
    % Validation of values
    if length(y) ~= size(X,1)
        error('The dimensions of y and X must agree.')
    end
    if mod(nmcmc,1) ~= 0 || nmcmc < 1
        error('nmcmc must be a postiive integer value.')
    end
    if mod(burn,1) ~= 0 || burn < 0
        error('burn must be a postive integer value.')
    end
    if mod(leafmin,1) ~= 0 || leafmin < 1
        error('leafmin must be a positive integer value.')
    end
    if gamma <= 0  || gamma >= 1
        error('gamma must be between 0 and 1')
    end
    if beta < 0
        error('beta must be positive.')
    end
    if p < 0 || p > 1
        error('p must be between 0 and 1')
    end
    if hottemp >= 1 || hottemp <= .01
        error('hottemp must be between .01 and 1')
    end
    if mod(swapfreq,1) ~= 0 || swapfreq < 1
        error('swapfreq must be an integer >= 1.')
    end
    % Add slash if necessary
    if ~strcmp(filepath(end),'/')
        filepath = strcat(filepath,'/');
    end
    % Create output directory or warn about overwriting
    if ~isdir(filepath)
        mkdir(filepath)
        disp(strcat(['NOTE: Creating output directory ',filepath]))
    else
        etext = strcat(['NOTE: May overwrite files in ',filepath]);
        disp(etext)
    end
    if suppress_errors_on_workers
        disp('NOTE: Errors suppressed on workers.');
    end
    
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
    swapaccepttotal = 0;
    swaptotal = 0;
    swapaccepttotal_global = 0;
    swaptotal_global = 0;
    
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    m = poolsize;
    % Harmonic Temperatures
    %delta = (1/hottemp - 1)/(m-1);
    %temps = [1, 1./(1 + delta*(1:m - 1))];
    
    % Sigmoidal temperatures
    j1 = log( 1/(-1 + 1.01) - 1);
    jm = log( 1/(-hottemp + 1.01) - 1);
    dm = (jm-j1)/(m-1);
    js = j1:dm:jm;
    temps = 1.01 - 1./(1 + exp(js));

    spmdsize = min([poolsize,m]);
    if spmdsize < 1
        error('Must have at least two processes to do parallel tempering.');
    end
       
    spmd(spmdsize)
        if parallelprofile
            mpiprofile on
        end
        if suppress_errors_on_workers
            oldwarnstate0 = warning('off','all');
        end
        % Suppress Matlab error for nearly singular matrix
        oldwarnstate = warning('off','MATLAB:nearlySingularMatrix');
        myname = labindex;
        master = 1; % master process labindex
        % Create independent Random Streams with a seed on each lab
        s = RandStream.create('mrg32k3a','Numstreams',m,...
            'StreamIndices',myname,'Seed',seed);
        RandStream.setGlobalStream(s);
        % Initialize root tree on each process
        mytemp = temps(myname);
        T = Tree(y,X,leafmin,gamma,beta,mytemp);
        if myname == master
            disp('Starting MCMC...')
        end
        for ii=1:(burn + nmcmc)
            % Propose a tree (with error handling)
            goodtree = 0;
            errcntr = 0;
            while ~goodtree
                try 
                    [Tstar,~,r,lr] = proposeTree(T,y,X,allprobs,p,mytemp,0);
                    goodtree = 1;
                catch ME % If error, try again
                    msg = getReport(ME);
                    warning(msg);
                end
                errcntr = errcntr + 1;
                if errcntr > 100
                    error('Possibly an infinite loop encountered.')
                end
            end
            if lr > log(rand)
                T = Tstar;
                naccept = naccept + 1;
                if r == 1
                    n_g_accept = n_g_accept + 1;
                    n_g_total = n_g_total + 1;
                    %tsize = tsize + 1;
                elseif r == 2
                    n_p_accept = n_p_accept + 1;
                    n_p_total = n_p_total + 1;
                    %tsize = tsize - 1;
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
            if mod(ii,swapfreq) == 0
                % Propose a switch of chains and send to all workers
                if myname == master
                    swapind = zeros(1,2);
                    swapind(1) = randsample(m-1,1);
                    swapind(2) = swapind(1) + 1;
                    swapind = labBroadcast(master,swapind);
                else
                    swapind = labBroadcast(master); 
                end
                % Send proposed swap to master
                if myname == swapind(1) && myname ~= master
                    labSend(T,master,1);
                    swaptotal = swaptotal + 1;
                end
                if myname == master
                    if myname ~= swapind(1)
                        Tstarswap1 = labReceive(swapind(1),1);
                    else
                        Tstarswap1 = T;
                        swaptotal = swaptotal + 1;
                    end
                end
                if myname == swapind(2) && myname ~= master
                    labSend(T,master,2);
                    swaptotal = swaptotal + 1;
                end
                if myname == master
                    if myname ~= swapind(2)
                        Tstarswap2 = labReceive(swapind(2),2);
                    else
                        Tstarswap2 = T;
                        swaptotal = swaptotal + 1;
                    end
                end
                if myname == master
                    swaptotal_global = swaptotal_global + 1;
                    lrswap = (Tstarswap2.Temp * Tstarswap1.Lliketree + Tstarswap1.Prior) + ...
                        (Tstarswap1.Temp * Tstarswap2.Lliketree + Tstarswap2.Prior) - ...
                        (Tstarswap1.Temp * Tstarswap1.Lliketree + Tstarswap1.Prior) - ...
                        (Tstarswap2.Temp * Tstarswap2.Lliketree + Tstarswap2.Prior);
                    if ~isfinite(lrswap)
                        error('Non-finite swap likelihood ratio.')
                    end
                    if lrswap > log(rand) % Accept
                        swapaccept = 1;
                    else
                        swapaccept = 0;
                    end
                    swapaccepttotal_global = swapaccepttotal_global + swapaccept;
                    swapaccept = labBroadcast(master,swapaccept);
                else
                    swapaccept = labBroadcast(master);
                end
                if swapaccept
                    if myname == master
                        if myname ~= swapind(1)
                            labSend(Tstarswap2,swapind(1))
                        else
                            T = Tstarswap2;
                            swapaccepttotal = swapaccepttotal + 1;
                        end
                        if myname ~= swapind(2)
                            labSend(Tstarswap1,swapind(2))
                        else
                            T = Tstarswap1;
                            swapaccepttotal = swapaccepttotal + 1;
                        end
                    elseif any(myname == swapind)
                        swapaccepttotal = swapaccepttotal + 1;
                        T = labReceive(master);
                    end
                    if any(myname == swapind) % Update temperature on Tree
                        T.Temp = mytemp;
                    end
                end
            end

            %if myname == master % Print progress
                if mod(ii,100) == 0
                    disp(['i = ',num2str(ii),', ID = ',num2str(myname),', llike = ',num2str(T.Lliketree),...
                        ', accept = ',num2str(naccept/ii),...
                        ', swapaccept = ',num2str(swapaccepttotal),'/',num2str(swaptotal),...
                        ', Size = ',num2str(T.Ntermnodes),...
                        ', temp = ',num2str(T.Temp)]);
                    if myname == master
                        disp(' ');
                    end
                end
            %end

            % Record Values       
            if ii > burn
                TREES{ii - burn} = thin_tree(T);
                treesize(ii - burn) = T.Ntermnodes;
                LLIKE(ii - burn) = T.Lliketree;
            end
            
        end
        perc_accept = naccept/(nmcmc + burn);
        move_accepts = [n_g_accept/n_g_total,...
            n_p_accept/n_p_total,...
            n_c_accept/n_c_total,...
            n_s_accept/n_s_total];
        swap_accept = swapaccepttotal/swaptotal;
        if k > 0
            C = unique(sort(LLIKE,'descend'),'rows','stable');
            if length(C) < k
                if ~saveall
                    if myname == master
                        % This note is only generated by the master.  Other
                        % workers will not display the message.
                        disp('NOTE: Number of unique trees is less than k. Less than k trees returned.')
                    end
                end
                k = length(C);
            end
            Treesub = cell(k,1);
            for mm=1:k
                I = find(C(mm) == LLIKE);
                Treesub{mm} = TREES{I(1)};
            end
            TREES = Treesub;
        end
        output = struct('Trees',{TREES},'llike',LLIKE,'acceptance',perc_accept,...
            'treesize',treesize,'move_accepts',move_accepts,...
            'swap_accept',swap_accept);
        % Save output
        if saveall
            savenames = 1:m;
        else
            savenames = master;
        end
        if any(myname == savenames)
            fname = strcat(filepath,'mcmc_id',num2str(myname),'.mat');
            swap_percent_global = swapaccepttotal_global/swaptotal_global;
            strt = tic;
            savedata(fname,output,swap_percent_global);
            stp = toc(strt);
            savetime = stp - strt;
        end
        % Turn on suppressed warnings
        warning(oldwarnstate);
        if suppress_errors_on_workers
            warning(oldwarnstate0);
            % warning('on','all')
        end
        if parallelprofile
            mpiprofile viewer
            mpiprofile off
        end
    end
end

function savedata(fname,output,swp_perc)
    save(fname,'output','swp_perc')
end

    
