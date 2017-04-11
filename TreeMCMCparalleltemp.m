function TreeMCMCparalleltemp(y,X,nmcmc,burn,leafmin,gamma,beta,p,hottemp,swapfreq,filepath)
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
    if spmdsize < m
        error('Number of temperatures cannot exceed number of processes.');
    end
    
    spmd(spmdsize)
        myname = labindex;
        master = 1; % master process labindex

        % Initialize root tree on each process
        mytemp = temps(myname);
        T = Tree(y,X,leafmin,gamma,beta,mytemp);
        for ii=1:(burn + nmcmc)
            [Tstar,~,r,lr] = proposeTree(T,y,X,allprobs,p,mytemp,0);
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
            
            %disp('here0!')
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
                %disp('here!')

                %if any(myname == swapind)
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
                %disp('here1')
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

    %             if myname == swapind(2)
    %                 labSend(T,master,2)
    %                 swaptotal = swaptotal + 1;
    %             elseif myname == master
    %                 Tstarswap2 = labReceive(swapind(2),2);
    %             end

                %disp('here2')

                if myname == master
                    swaptotal_global = swaptotal_global + 1;
%                     temp1 = Tstarswap1.Temp;
%                     llike1 = Tstarswap1.Lliketree;
%                     temp2 = Tstarswap2.Temp;
%                     llike2 = Tstarswap2.Lliketree;
                    lrswap = (Tstarswap2.Temp * Tstarswap1.Lliketree + Tstarswap1.Prior) + ...
                        (Tstarswap1.Temp * Tstarswap2.Lliketree + Tstarswap2.Prior) - ...
                        (Tstarswap1.Temp * Tstarswap1.Lliketree + Tstarswap1.Prior) - ...
                        (Tstarswap2.Temp * Tstarswap2.Lliketree + Tstarswap2.Prior);
                    %[temp1,llike1,temp2,llike2,lrswap]
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

                %disp('here3')

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
                        ', Size = ',num2str(T.Ntermnodes)]);
                    if myname == master
                        disp(' ');
                    end
                end
            %end
            
            % Record Values       
            if ii > burn
                TREES{ii - burn} = T;
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

        output = struct('Trees',{TREES},'llike',LLIKE,'acceptance',perc_accept,...
            'treesize',treesize,'move_accepts',move_accepts,...
            'swap_accept',swap_accept);
        
        fname = strcat(filepath,'mcmc_id',num2str(myname),'.mat');
        swap_percent_global = swapaccepttotal_global/swaptotal_global;
        strt = tic;
        savedata(fname,output,swap_percent_global);
        stp = toc(strt);
        savetime = stp - strt
        
        % Save Output
        
        % Keep only the true chain
        % output = output{1};
    end
    %swap_percent_global = swapaccepttotal_global{1}/swaptotal_global{1};
    %output{'swap_accept_global'} = swapaccepttotal_global{1}/swaptotal_global{1};
    %output = output0{1};
end

function savedata(fname,output,swp_perc)
    save(fname,'output','swp_perc')
end

    