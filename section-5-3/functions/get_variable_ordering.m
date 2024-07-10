% Get the ordered list of variable names from their IDs.
%
% Notes: 
%   - The ordering is based on the decision rule ordering 
%     in the Dynare model.
%
% Inputs:
%   M_       - Returned by Dynare
%   options_ - Returned by Dynare
%   oo_      - Returned by Dynare
%   vars     - Vector of variable IDs for which the names are to be retrieved.
%
% Outputs:
%   ordering - Column vector of ordered variable names corresponding to the input IDs.
%
function ordering=get_variable_ordering(M_, options_, oo_, vars)
  cellordering = M_.endo_names(oo_.dr.order_var(vars));
  ordering = string(cell2mat(cellordering(1)));
  for i=2:length(cellordering)
    ordering = [ordering; string(cell2mat(cellordering(i)))];
  end
end
