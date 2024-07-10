% Compute the transmission effect of a channel that goes through x in some period.
%
% Inputs:
%   M_      - Returned by Dynare
%   B       - Returned by `varma_to_static`. Corrsponds to B in the static 
%             representation of Wegner etal (2024).
%   Oomega  - Returned by `varma_to_static`. Corrsponds to Omega in the static 
%             representation of Wegner etal (2024).
%   idx_x   - Index of variable through which the channel should go in some period.
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
function transmission_effect=through_x_some_period(M_, B, Oomega, idx_x, k)

  B_tilde = B;
  Oomega_tilde = Oomega;

  % Working through negation. First find effect never going through x
  % and then negate to find effect going through x in some period. 
  B_tilde(idx_x:k:end, :) = 0;
  Oomega_tilde(idx_x:k:end, :) = 0;

  irfs_never = irf_static_model(M_, B_tilde, Oomega_tilde, k);
  irfs = irf_static_model(M_, B, Oomega, k);

  transmission_effect = irfs - irfs_never; 
end



