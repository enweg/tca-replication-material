% Compute the transmission effect of a channel that can go through one 
% variable (x) but not through y or z.
%
% Inputs:
%   M_      - Returned by Dynare
%   B       - Returned by `varma_to_static`. Corrsponds to B in the static 
%             representation of Wegner etal (2024).
%   Oomega  - Returned by `varma_to_static`. Corrsponds to Omega in the static 
%             representation of Wegner etal (2024).
%   idx_x   - Variable index for the variable through which the channel can go.
%   idx_y   - Variable index for the variable through which the channel cannot go.
%   idx_z   - Variable index for the variable through which the channel cannot go.
%   k       - Integer specifying the number of variables in the system.
%   h       - Vector of horizons.
%
% Outputs:
%   transmission_effect - 3D array of transmission effects (k x m x horizon+1), 
%                         where horizon is determined based on the size of B and k.
%
% References: 
%   - Wegner, E., Lieb, L., Smeekes, S., & Wilms, I. (2024). 
%     Transmission Channel Analysis in Dynamic Models. 
%     arXiv preprint arXiv:2405.18987.
function transmission_effect=through_x_not_y_not_z(M_, B, Oomega, idx_x, idx_y, idx_z, k, h)

  if length(h) > 1
    transmission_effect = through_x_not_y_not_z(M_, B, Oomega, idx_x, idx_y, idx_z, k, h(1));
    for i=2:length(h)
      tmp = through_x_not_y_not_z(M_, B, Oomega, idx_x, idx_y, idx_z, k, h(i));
      transmission_effect = cat(3, transmission_effect, tmp);
    end
    return
  end

  % First term: Q[and_{t=0}^h not y_t and not z_t]
  B_tilde = B; 
  Oomega_tilde = Oomega;
  for t=0:h
    B_tilde(idx_y+(t*k), :) = 0; 
    Oomega_tilde(idx_y+(t*k), :) = 0;
    B_tilde(idx_z+(t*k), :) = 0; 
    Oomega_tilde(idx_z+(t*k), :) = 0;
  end
  term_1 = irf_static_model(M_, B_tilde, Oomega_tilde, k);
  
  % Second term: Q[and_{t=0}^h not x_t and not y_t and not z_t]
  B_tilde = B; 
  Oomega_tilde = Oomega;
  for t=0:h
    B_tilde(idx_y+(t*k), :) = 0; 
    Oomega_tilde(idx_y+(t*k), :) = 0;
    B_tilde(idx_z+(t*k), :) = 0; 
    Oomega_tilde(idx_z+(t*k), :) = 0;
    B_tilde(idx_x+(t*k), :) = 0;
    Oomega_tilde(idx_x+(t*k), :) = 0;
  end
  term_2 = irf_static_model(M_, B_tilde, Oomega_tilde, k);

  transmission_effect = term_1 - term_2;
  transmission_effect = transmission_effect(:, :, h+1);
end



