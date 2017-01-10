function out = randsub(S)
    % Randomly partitions splits the groups uniformly 
    % Gives the left tree rule.
    % S must be a unique one-columned table
    S = table2cell(S);
    ns = size(S,1);
    nu = floor(ns/2);
    groupsizes = zeros(nu,1);
    for ii = 1:nu
        groupsizes(ii) = nchoosek(ns,ii);
    end
    cumgroup = cumsum(groupsizes);
    gind = randsample(sum(groupsizes),1);
    numgrpL = find(gind <= cumgroup,1);
    grps = nchoosek(S,numgrpL);
    out = grps(randsample(size(grps,1),1),:);
end

