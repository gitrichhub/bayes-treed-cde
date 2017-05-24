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


function m = mfunc(Tind,Tmultiset)
    m = 1;
    % If the likelihood matches, they may be the same tree.
    llike = Tmultiset{Tind}.Lliketree;
    for ii = find(1:length(Tmultiset) ~= Tind)
        Ttemp = Tmultiset{Tind};
        if Ttemp.Lliketree == llike
             % Check Tree size
            if length(Ttemp.Allnodes) == length(Tmultiset{ii}.Allnodes)
                tm = treematch(Ttemp,Tmultiset{ii},[],[]);
                m = m + tm;
            end
        end
    end
end