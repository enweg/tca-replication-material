% Convert a cellarray used by Dynare into a standard vector
%
% Inputs:
%   ca - Cell array where each cell contains a string. This is typically the format 
%        used by Dynare for storing names or labels.
%
% Outputs:
%   v  - Vector of strings, converted from the input cell array.
function v=dynare_cellarray_to_vec(ca)
  v = repelem("", length(ca));
  for i=1:length(ca)
    v(i) = string(cell2mat(ca(i)));
  end
  v = vec(v);
end
