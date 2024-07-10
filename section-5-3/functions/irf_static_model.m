% Compute the Impulse Response Functions (IRFs) in the static representation 
% x = B*x + Oomega*eps of Wegner etal (2024)
%
% Inputs:
%   M_      - Returned by Dynare
%   B       - See Wegner etal (2024). Obtained using `varma_to_static`
%   Oomega  - See Wegner etal (2024). Obtained using `varma_to_static`
%   k       - Number of variables in the system. Corresponds to the 
%             dimension of A0 returned by `get_varma_coeffs`
%
% Outputs:
%   Ttheta  - 3D array of IRFs. The dimensions are (k, nexo, horizon+1), 
%   where k is the number of variables, nexo is the number of exogenous shocks, 
%   and horizon+1 is the number of periods for the IRFs. In general, nexo=k.
%
% References: 
%   - Wegner, E., Lieb, L., Smeekes, S., & Wilms, I. (2024). 
%     Transmission Channel Analysis in Dynamic Models. 
%     arXiv preprint arXiv:2405.18987.
function Ttheta=irf_static_model(M_, B, Oomega, k)
  % k is the number of variables in the system including observed and
  % state variables. 
  nexo = M_.exo_nbr;
  Ttheta_static = inv(eye(size(B, 1)) - B) * Oomega; 
  horizon = floor(size(B, 1) / k - 1);
  Ttheta = zeros(k, nexo, horizon+1);
  for h=0:horizon
    Ttheta(:, :, h+1) = Ttheta_static((h*k+1):((h+1)*k), 1:nexo);
  end
end

