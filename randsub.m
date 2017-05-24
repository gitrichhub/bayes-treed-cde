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

