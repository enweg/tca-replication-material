% Get Dynare variable IDs for given variable names.
%
% Inputs:
%   M_        - Returned by Dynare
%   options_  - Returned by Dynare
%   oo_       - Returned by Dynare
%   varnames  - Vector of strings, where each string is a variable name
%               as specified in the .mod file.
%
% Outputs:
%   obs_var   - Vector of IDs corresponding to the given variable names. 
%               These IDs are the internal Dynare indices for the specified 
%               variables.
%
% Notes:
%   - The function assumes that all variable names provided in varnames are 
%     present in the model's endogenous variable names.
%   - The IDs are determined based on the position of each variable name 
%     within the endogenous variable names list (M_.endo_names) and then mapped 
%     to the order used by Dynare's decision rules (oo_.dr.inv_order_var).
function obs_var=get_obsvar_ids(M_, options_, oo_, varnames)

  ids = repelem(1, length(varnames));
  endoname = dynare_cellarray_to_vec(M_.endo_names);
  for i=1:length(varnames)
    name = varnames(i);
    ids(i) = find(endoname == name);
  end
     
  obs_var = oo_.dr.inv_order_var(ids);
end
