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

function [out,llikes,priors] = getpost(Tmultiset,X)
    K = length(Tmultiset);
    priors = zeros(K,1);
    llikes = zeros(K,1);
    for ii=1:K
        llikes(ii) = Tmultiset{ii}.Lliketree;
        % TODO: Make this more efficient with an indicator variable!
        priors(ii) = prior_eval(Tmultiset{ii},X);
    end
    post = llikes + priors;
    a = max(post);
    out = a + log(sum(exp(post - a)));
end