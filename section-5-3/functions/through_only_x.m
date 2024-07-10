% Compute the transmission effect of a channel that can go only through 
% a single variable at some horizon. 
%
% Inputs:
%   M_      - Returned by Dynare
%   B       - Returned by `varma_to_static`. Corrsponds to B in the static 
%             representation of Wegner etal (2024).
%   Oomega  - Returned by `varma_to_static`. Corrsponds to Omega in the static 
%             representation of Wegner etal (2024).
%   idx_x   - Variable index for the variable through which the channel can go.
%   k       - Integer specifying the number of variables in the system.
%
% Outputs:
%   transmission_effect - 3D array of transmission effects (k x m x horizon+1), 
%                         where horizon is determined based on the size of B and k.
%
% References: 
%   - Wegner, E., Lieb, L., Smeekes, S., & Wilms, I. (2024). 
%     Transmission Channel Analysis in Dynamic Models. 
%     arXiv preprint arXiv:2405.18987.
function transmission_effect=through_only_x(M_, B, Oomega, idx_x, k)
  not_through_idx = setdiff(1:k, idx_x);
  transmission_effect = through_not_x(M_, B, Oomega, not_through_idx, k);
end
