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

classdef SplitRule
    properties
        Varcol
        Varclass
        Varrule        
    end
    methods
        % WE DON'T NEED TO FEED IN THE ENTIRE X! CHANGE FOR LATER!
        function obj = SplitRule(X,colnum,rule)
            if nargin == 3
                if isa(X,'table')
                    obj.Varcol = colnum;
                    vartypes = varfun(@class,X,'OutputFormat','cell');
                    obj.Varclass = vartypes(colnum);
                    if strcmp(obj.Varclass,'double')
                        if isa(rule,'double')
                            obj.Varrule = rule;
                        else
                            error('Rule must match be a double')
                        end
                    elseif strcmp(obj.Varclass,'cell')
                        if isa(rule,'cell')
                            obj.Varrule = rule;
                        else
                            error('Rule must be a cell for factors')
                        end
                    else
                        error('Data columns must be either doubles or cell')
                    end
                else
                    error('Design matrix must be of class "table"');
                end
            elseif nargin ~= 0 
                error('Constructor takes 3 arguments')
            end
        end
        % Add in functions for changing rules, columns, etc.
    end
    methods (Static)
        % Determine if two rules are identical
        function out = equalrule(rule1,rule2)
            if ~isempty(rule1) %~isempty(rule1.Varcol) && ~isempty(rule1.Varrule)
                if ~isempty(rule2) %~isempty(rule2.Varcol) && ~isempty(rule2.Varrule)
                    if (rule1.Varcol == rule2.Varcol) && (rule1.Varrule == rule2.Varrule)
                        out = 1;
                    else
                        out = 0;
                    end
                else % rule 1 is not empty, rule 2 is empty
                    out = 0;
                end
            else % rule1 is empty
                if ~isempty(rule2) %~isempty(rule2.Varcol) && ~isempty(rule2.Varrule) % rule2 is not empty
                    out = 0;
                elseif isempty(rule2) %isempty(rule2.Varcol) && isempty(rule2.Varrule)
                    out = 1; % both are empty;
                end
            end     
        end
    end
end