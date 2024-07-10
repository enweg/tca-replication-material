% Get VARMA coefficients for the VARMA representation of a DSGE model.
%
% Inputs:
%   M_        - Returned by Dynare
%   options_  - Returned by Dynare
%   oo_       - Returned by Dynare
%   p         - The VAR order of the VARMA model determined using 
%               `determine_varma_order`.
%   q         - The MA order of the VARMA model determined using 
%               `determine_varma_order`.
%   obs_var   - (Optional) Vector of observed variables. If not provided, the function 
%               will derive this from the given DSGE model structures.
%
% Outputs:
%   A0        - Matrix representing the contemporaneous impact matrix in the VARMA model.
%   Phis      - Array containing the VAR coefficients of the VARMA model.
%   Psis      - Array containing the MA coefficients of the VARMA model.
%   vars      - Vector of variable IDs corresponding to the VARMA representation.
%
% Notes:
%   - If p <= 2 and q <= 1, the function uses Corollary 1 or 2 from the referenced paper
%     to calculate the coefficients directly.
%   - For general cases, the function uses Proposition 1, which requires that p - 1 = q.
%     An error is raised if this condition is not satisfied or if D is not invertible.
%
% References: 
%   - Morris, S. D. (2016). VARMA representation of DSGE models. 
%     Economics Letters, 138, 30–33. https://doi.org/10.1016/j.econlet.2015.11.027
function [A0, Phis, Psis, vars]=get_varma_coeffs(M_, options_, oo_, p, q, obs_var)

  if nargin==5
    [Ax, Ay, Bx, Cx, Cy, D, vars] = get_ABCD_varma(M_, options_, oo_);
  elseif nargin==6
    [Ax, Ay, Bx, Cx, Cy, D, vars] = get_ABCD_varma(M_, options_, oo_, obs_var);
  end


  n = size(Cx, 1);
  m = size(Ax, 1);
  Phis = zeros(n, n, p);
  Psis = zeros(n, n, q);

  vars = vars((end - n + 1):end);
  vars = get_variable_ordering(M_, options_, oo_, vars);

  % Main assumption
  if rank(D) < size(D, 1)
    error("D is not invertible.");
  end

  if p<=2 && q<=1
    % case of corollary 1 or 2
    % double checking rank of Cx
    if rank(Cx) ~= size(Cx, 1)
      error("Cannot have p<=2 and q<=1 if Cx is not full rank.")
    end
    Phis(:, :, 1) = D \ ((Cx * Ax) / Cx + Cy);
    if p==2
      Phis(:, :, 2) = D \ (Ay - ((Cx * Ax) / Cx) * Cy);
    end
    Psis(:, :, 1) = D \ ((Cx * Bx / D - Cx * Ax / Cx)*D);
    A0 = D \ eye(n);
    return;
  end

  % case of general proposition 1
  if p - 1 ~= q
    error("Proposition states that p-1=q. This is not satisfied.");
  end
  kappa = p-2;

  C = [Cx Cy];
  A = [Ax Ay; Cx Cy];
  B = [Bx; D];

  F = C; 
  for k=1:kappa
    F = [C*A^k; F];
  end
  Fx = F(:, 1:m);
  Fx_plus = pinv(Fx);
  Fy = F(:, (m+1):end);
  
  Gs = zeros(n*(kappa+1), n, kappa+1);
  Gi = D;
  Gs(:, :, 1) = [Gi; zeros(n*kappa, n)];
  Gi = [C*B; Gi];
  Gs(:, :, 2) = [Gi; zeros(n*(kappa-1), n)];
  for i=3:(kappa+1)
    Gi = [C*A^(i-2)*B; Gi];
    Gs(:, :, i) = [Gi, zeros(n*(kappa+1-i), n)];
  end
    
  coef_et = (Fx_plus * Gs(:, :, 1));
  pinv_coef_et = pinv(coef_et);
  A0 =  pinv_coef_et * Fx_plus(:, 1:n); 
  Phis(:, :, end) = pinv_coef_et *  (Ay - Ax * Fx_plus * Fy);
  Phis(:, :, end-1) = pinv_coef_et * (Ax * Fx_plus(:, (end-n+1):end) + Fx_plus * Fy);
  for k=1:kappa
    Phis(:, :, kappa) = pinv_coef_et * (Ax * Fx_plus(:, ((k-1)*n+1):(k*n)) - Fx_plus(:, (k*n+1):((k+1)*n)));
  end
  
  Psis(:, :, 1) = pinv_coef_et * (Bx - Ax * Fx_plus*Gs(:, :, end));
  for k=1:kappa
    Psis(:, :, k) = pinv_coef_et * (Fx_plus * Gs(:, :, k+1) - Ax * Fx_plus * Gs(:, :, k));
  end


end


