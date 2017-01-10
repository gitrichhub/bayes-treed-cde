function TreeMCMC(y,X,nmcmc,burn)
    % Check for format of input variables
    if ~isnumeric(y)
        error('y must be numeric.')
    end
    if ~istable(X)
        error('X must be a table.')
    end
    if length(nmcmc) ~= 1 || ~isinteger(nmcmc)
        error('input arg "nmcmc" must be a scalar integer.')
    end
    if length(burn) ~= 1 || ~isinteger(burn)
        error('input arg "burn" must be a scalar integer.')
    end
    
    
    % Initialize root tree
    
    
    
    % Posterior trees
    TREES = cell(nmcmc,1);
    for ii=1:(burn + nmcmc)
        % Initialize 
        
        
        while ii > burn
           
           
        end
    end






end


