% Compute the transmission effect of a channel not going through any 
% variables in `idx_x` for all horizons.
%
% Inputs:
%   M_      - Returned by Dynare
%   B       - Returned by `varma_to_static`. Corrsponds to B in the static 
%             representation of Wegner etal (2024).
%   Oomega  - Returned by `varma_to_static`. Corrsponds to Omega in the static 
%             representation of Wegner etal (2024).
%   idx_x   - Vector of indices specifying the variables to exclude. 
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
function transmission_effect=through_not_x(M_, B, Oomega, idx_x, k)
  B_tilde = B;
  Oomega_tilde = Oomega; 

  for i=1:length(idx_x)
    ix = idx_x(i);
    B_tilde(ix:k:end, :) = 0;
    Oomega_tilde(ix:k:end, :) = 0;
  end

  transmission_effect = irf_static_model(M_, B_tilde, Oomega_tilde, k);
end

